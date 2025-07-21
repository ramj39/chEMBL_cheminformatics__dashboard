import csv
from datetime import datetime

def collect_feedback(user_input: str) -> dict:
    """Validate and clean feedback input."""
    feedback = user_input.strip()
    return {"feedback": feedback, "status": "submitted" if feedback else "empty"}

#def save_feedback(feedback_data: dict, filepath="feedback_log.csv") -> None:
#    """Save feedback to a local CSV file with timestamp."""
#   timestamp = datetime.now().isoformat(timespec="seconds")
#   with open(filepath, "a", newline="", encoding="utf-8") as f:
#       writer = csv.writer(f)
#       writer.writerow([timestamp, feedback_data["status"], feedback_data["feedback"]])
import csv
from datetime import datetime

def save_feedback(feedback_data: dict, filepath="feedback_log.csv") -> None:
    timestamp = datetime.now().isoformat(timespec="seconds")
    status = feedback_data.get("status", "submitted")
    comment = feedback_data.get("feedback", "")

    with open(filepath, "a", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([timestamp, status, comment])

#def save_feedback(feedback_data: dict, filepath="feedback_log.csv") -> None:
#   from datetime import datetime
#   import csv

#   timestamp = datetime.now().isoformat(timespec="seconds")
#   status = feedback_data.get("status", "submitted")
#   comment = feedback_data.get("feedback", "")

#   with open(filepath, "a", newline="", encoding="utf-8") as f:
#       writer = csv.writer(f)
#       writer.writerow([timestamp, status, comment])

