"""Langgraph-based MaizeChatBot."""
import json
from operator import add
from typing import List, Dict, Any, Optional
from typing_extensions import TypedDict, Annotated
from langgraph.graph import StateGraph, START, END
from langgraph.graph.message import add_messages
from langgraph.checkpoint.memory import MemorySaver

from utils import (
    routing_client, neo4j_driver, hybrid_retrieval,
    cfg, logger, language_adapter, generate_answer,
)

class GraphState(TypedDict):
    """Graph state."""
    question: str
    session_id: str
    
    task: Optional[str]
    entity: Optional[str]
    rewritten_question: Optional[str]
    history_summary: Optional[str]
    
    history: Annotated[List[Dict[str, str]], add]
    turn_count: int
    last_entity: Optional[str]
    last_task: Optional[str]

    result: Optional[Dict[str, Any]]

async def summarize_node(state: GraphState):
    """This node summarizes the recent conversation.
    
    Args:
        state (GraphState): The graph state.
    """
    history = state.get("history", [])
    if not history:
        return {"history_summary": "No prior conversation."}

    recent = []
    for msg in history[-4:]:
        role = "User" if msg["role"] == "user" else "Assistant"
        content = msg["content"][:150]
        if "The information about" in content:
            content = f"[KG Result: {content.replace('The information about ', '')}]"
        recent.append(f"{role}: {content[:100]}")
    dialogue = " | ".join(recent)
    prompt = cfg['prompts']['summary_prompt'].format(dialogue=dialogue)

    try:
        resp = routing_client.chat.completions.create(
            model=cfg['routing_model']['model'],
            messages=[{"role": "user", "content": prompt}],
            temperature=0.1, max_tokens=64,
        )
        return {"history_summary": resp.choices[0].message.content.strip()}
    except Exception as e:
        logger.warning(f"Summarize failed: {e}")
        return {"history_summary": f"Most recently discussed: {state.get('last_entity', 'unknown')}"}
    

async def rewrite_node(state: GraphState):
    prompt = cfg['prompts']['rewrite_prompt'].format(
        history_summary=state.get("history_summary", ""),
        question=state["question"],
    )
    try:
        resp = routing_client.chat.completions.create(
            model=cfg['routing_model']['model'],
            messages=[{"role": "user", "content": prompt}],
            temperature=0.1, max_tokens=80,
        )
        return {"rewritten_question": resp.choices[0].message.content.strip().strip('"\': \n')}
    except Exception:
        return {"rewritten_question": state["question"]}
    
async def route_node(state: GraphState):
    q = state.get("rewritten_question") or state["question"]
    try:
        resp = routing_client.chat.completions.create(
            model=cfg['routing_model']['model'],
            messages=[
                {"role": "system", "content": cfg['prompts']['task_judgement']},
                {"role": "user", "content": q},
            ],
            temperature=0.01, max_tokens=128,
        )
        content = resp.choices[0].message.content.strip()
        if "```json" in content:
            content = content.split("```json")[1].split("```")[0].strip()
        decision = json.loads(content)
    except Exception as e:
        logger.error(f"Routing failed: {e}")
        decision = {"task": "RAG", "entity": None}

    return {
        "task": decision.get("task", "RAG"),
        "entity": decision.get("entity"),
        "rewritten_question": q,
    }

async def kg_node(state: GraphState):
    variety = language_adapter.normalize_entity(state["entity"])
    CYPHER = cfg['cypher']['query_cypher']
    try:
        with neo4j_driver.session() as sess:
            record = sess.run(CYPHER, name=variety).single()
            if not record or not record.get("all_meta"):
                result = {"type": "KG", "status": "not_found",
                         "result": {"meta_dict": {}, "triples": []}}
            else:
                meta_dict = {item["name"]: item for item in record["all_meta"] if item.get("name")}
                result = {
                    "type": "KG", "status": "ok",
                    "result": {"meta_dict": meta_dict, "triples": record["triples"]},
                }
    except Exception as e:
        logger.error(f"KG error: {e}")
        result = {"type": "KG", "status": "error",
                  "result": {"meta_dict": {}, "triples": []}}

    return {
        "result": result,
        "last_entity": variety,
        "last_task": "KG",
        "turn_count": state["turn_count"] + 1,
        "history": [
            {"role": "user", "content": state["question"]},
            {"role": "assistant", "content": f"The information about {variety}"},
        ],
    }

async def rag_node(state: GraphState):
    question = state["rewritten_question"] or state["question"]
    lang = language_adapter.detect_language(state["question"])
    rag_query = language_adapter.normalize_rag_query(question, lang)

    retrieved = hybrid_retrieval(rag_query)
    sources = [{"source_id": h["source_id"]} for h in retrieved]
    context_text = "\n".join(
        [f"[{i+1}] {h['text'][:200]}..." for i, h in enumerate(retrieved[:3])]
    )

    messages = [
        {"role": "system", "content": cfg['prompts']['rag_prompt']},
        *list(state.get("history", [])),
        {"role": "user", "content": f"Relevant excerpts:\n{context_text}\n\nQuestion: {rag_query}"},
    ]
    raw = generate_answer(messages, temperature=0.3)
    answer = raw.split("</think>")[-1].strip() if "</think>" in raw else raw.strip()
    answer = language_adapter.localize_rag_answer(answer, lang)

    return {
        "result": {"type": "RAG", "answer": answer, "sources": sources},
        "last_task": "RAG",
        "turn_count": state["turn_count"] + 1,
        "history": [
            {"role": "user", "content": state["question"]},
            {"role": "assistant", "content": answer},
        ],
    }

def need_context_branch(state: GraphState):
    return "summarize" if state.get("turn_count", 0) > 0 else "route"

def task_branch(state: GraphState):
    return "kg" if (state["task"] == "KG" and state.get("entity")) else "rag"

def _build_graph():
    g = StateGraph(GraphState)

    g.add_node("summarize", summarize_node)
    g.add_node("rewrite", rewrite_node)
    g.add_node("route", route_node)
    g.add_node("kg", kg_node)
    g.add_node("rag", rag_node)

    g.add_conditional_edges(START, need_context_branch,
                            {"summarize": "summarize", "route": "route"})
    g.add_edge("summarize", "rewrite")
    g.add_edge("rewrite", "route")
    g.add_conditional_edges("route", task_branch, {"kg": "kg", "rag": "rag"})
    g.add_edge("kg", END)
    g.add_edge("rag", END)

    return g.compile(checkpointer=MemorySaver())

graph_app = _build_graph()
