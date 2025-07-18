# feedback.py

def collect_feedback(user_input: str) -> dict:
    """Validate and clean feedback input."""
    feedback = user_input.strip()
    return {"feedback": feedback, "status": "submitted" if feedback else "empty"}

def save_feedback(feedback_data: dict, filepath="feedback_log.txt") -> None:
    """Save feedback to a local file."""
    with open(filepath, "a", encoding="utf-8") as f:
        f.write(feedback_data["feedback"] + "\n")
