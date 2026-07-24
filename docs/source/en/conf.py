# English documentation tree -- see ../_shared/conf_base.py for shared settings.
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_shared"))

from conf_base import *  # noqa: F401,F403,E402

language = "en"
html_title = f"sqc {release}"  # noqa: F405

# Language switcher target (consumed by _templates/sidebar/language-switch.html)
html_context = {
    "current_language": "en",
    "other_language_name": "中文",
    "other_language_url": "../zh/",
}
