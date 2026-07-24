# 中文文档树 -- 共享配置见 ../_shared/conf_base.py。
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "_shared"))

from conf_base import *  # noqa: F401,F403,E402

language = "zh_CN"
html_title = f"sqc {release}"  # noqa: F405

# 语言切换目标(供 _templates/sidebar/language-switch.html 使用)
html_context = {
    "current_language": "zh",
    "other_language_name": "English",
    "other_language_url": "../",
}
