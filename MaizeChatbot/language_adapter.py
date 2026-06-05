import re
import yaml
from openai import OpenAI

def contains_chinese(text: str) -> bool:
    """
    Check if string contains chinese
    """
    return bool(re.search(f'[\u4e00-\u9fff]', text))

class LanguageAdapter:
    def __init__(self, cfg: dict):
        self.cfg = cfg
        self.client = OpenAI(
            api_key=cfg['routing_model']['api_key'],
            base_url=cfg['routing_model']['base_url']
        )
        
        # sample alias
        self.variety_alias = self._load_variety_alias()
        
    def _load_variety_alias(self) -> dict:
        """
        Load en->zh variety alias table
        """
        path = self.cfg.get("variety_alias_path")
        if not path:
            return {}
        
        try:
            with open(path, 'r', encoding="utf-8") as f:
                return yaml.safe_load(f) or {}
        except Exception:
            return {}
        
    # language detection
    def detect_language(self, text: str) -> str:
        return "zh" if contains_chinese(text) else "en"
    
    def translate(self, text: str, src: str, tgt: str) -> str:
        """
        Generic translator using cloud LLM
        """
        if src == tgt:
            return text
        
        prompt = (
            f"You are a professional scientific translator specializing in maize genetics "
            f"and plant breeding research. Translate the following text strictly "
            f"from {src} to {tgt}.\n\n"
            f"Requirements:\n"
            f"1. Preserve the original question structure and intent without reinterpretation.\n"
            f"2. Maintain the original level of abstraction (e.g., breeding pipeline, "
            f"methodology, genetic mechanism, or biological function).\n"
            f"3. Use terminology standard in maize genetics and doubled haploid (DH) breeding literature.\n"
            f"4. Do NOT add explanations, examples, or change the scientific perspective.\n\n"
            f"Text:\n{text}"
        )

        
        resp = self.client.chat.completions.create(
            model=self.cfg['routing_model']['model'],
            messages=[{"role": "user", "content": prompt}],
            temperature=0.0,
            max_tokens=512
        )
        return resp.choices[0].message.content.strip()
    
    # RAG helper
    def normalize_rag_query(self, question: str, lang: str) -> str:
        """
        Normalize RAG input to english 
        """
        if lang == "zh":
            return self.translate(question, src="Chinese", tgt="English")
        return question
    
    def localize_rag_answer(self, answer: str, lang: str) -> str:
        """
        Localize RAG answer to user's language
        """
        if lang == "zh":
            return self.translate(answer, src="English", tgt="Chinese")
        return answer
    
    # KG helper
    def normalize_entity(self, entity: str) -> str:
        """
        Normalize KG entity name (EN -> ZH if alias exists)
        """
        if not entity:
            return entity
        return self.variety_alias.get(entity, entity)
