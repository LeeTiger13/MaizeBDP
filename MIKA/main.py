import uuid
from fastapi import FastAPI, Request
from fastapi.responses import StreamingResponse
from fastapi.middleware.cors import CORSMiddleware
from utils import logger
from mika_langgraph import graph_app


app = FastAPI(title="MaizeBDP Hybrid Gateway")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])

@app.post("/chat")
async def chat(request: Request):
    data = await request.json()
    question = data.get("question")
    session_id = data.get("session_id") or str(uuid.uuid4())

    config = {"configurable": {"thread_id": session_id}}
    snapshot = graph_app.get_state(config).values or {}

    inputs = {
        "question": question,
        "session_id": session_id,
        "turn_count": snapshot.get("turn_count", 0),
        "history": snapshot.get("history", []),
        "last_entity": snapshot.get("last_entity"),
        "last_task": snapshot.get("last_task"),
    }

    final_state = await graph_app.ainvoke(inputs, config=config)
    result = final_state["result"]
    result["session_id"] = session_id
    logger.info(f"[{session_id}] task={final_state.get('task')} entity={final_state.get('entity')}")
    return result

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8082)
