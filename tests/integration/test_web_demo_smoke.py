"""tests.integration.test_web_demo_smoke — web_demo_v2 build smoke test.

Ensures the sqc-based Gradio demo imports and constructs its UI without
launching a server, and that experimental tabs stay hidden in v1 (D3).
Skipped when gradio (a demo-only optional dependency) is not installed.
"""
import pytest

pytest.importorskip("gradio", reason="gradio is a demo-only optional dependency")


def test_web_demo_v2_builds_without_launch():
    import web_demo_v2 as demo

    # v1: experimental tabs (e.g. Z-Crosstalk) are hidden by default.
    assert demo.SHOW_EXPERIMENTAL is False

    app = demo.build_app()
    # Gradio Blocks object constructed successfully (no .launch()).
    assert app.__class__.__name__ == "Blocks"
