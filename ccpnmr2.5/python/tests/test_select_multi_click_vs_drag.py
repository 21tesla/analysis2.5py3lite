"""
Regression tests for WindowFrame.selectMulti's click-vs-drag classification.

x0/x1/y0/y1 arrive as 0..1 *fractions* of the canvas, never as pixels.
Commit f2360234 added `abs(fraction - fraction) < 5` as a "click tolerance",
but the delta of two fractions is always < 1, so the test was always true and
*every* drag was routed to selectSingle -- examineRegion (which runs
searchPeaks/findPeaks) was never reached, so box-select and box-pick silently
did nothing. The fix scales the fraction deltas back to pixels before
applying the 5px tolerance.

These tests pin the boundary so the unit mismatch cannot return.
"""

from unittest.mock import MagicMock

from ccpnmr.analysis.frames.WindowFrame import WindowFrame


def _make_canvas(width=1000, height=600):
    canvas = MagicMock()
    canvas.winfo_width.return_value = width
    canvas.winfo_height.return_value = height
    canvas.handler = None
    return canvas


def _make_self_mock():
    self_mock = MagicMock()
    self_mock.topPopup = MagicMock()
    self_mock.unpostMenu = MagicMock()
    self_mock.selectSingle = MagicMock()
    self_mock.examineRegion = MagicMock()
    self_mock.setCurrentObjects = MagicMock()
    # canvas.prev_cross / prev_box must not be "present" so the hasattr guards
    # in selectMulti do not try to del attributes on a plain MagicMock.
    self_mock.canvas_has = False
    return self_mock


def _call(self_mock, canvas, x0, y0, x1, y1, state=0):
    WindowFrame.selectMulti(
        self_mock,
        canvas,
        a0=0.1, b0=0.2, a1=0.3, b1=0.4,
        x0=x0, y0=y0, x1=x1, y1=y1,
        state=state,
    )


def test_treated_as_click_when_drag_smaller_than_5px():
    # 3px x 2px drag on a 1000x600 canvas -> click.
    self_mock = _make_self_mock()
    canvas = _make_canvas(1000, 600)
    _call(self_mock, canvas, x0=0.100, y0=0.500, x1=0.103, y1=0.503)  # 3px, 1.8px
    self_mock.selectSingle.assert_called_once()
    self_mock.examineRegion.assert_not_called()


def test_treated_as_drag_when_x_moves_far_in_pixels_but_small_in_fraction():
    # The bug: a big 3.5-unit px drag is a tiny fraction (0.0035).
    self_mock = _make_self_mock()
    canvas = _make_canvas(1000, 600)
    # 2px x 200px drag: y moves 200px.
    _call(self_mock, canvas, x0=0.100, y0=0.500, x1=0.102, y1=0.833)
    self_mock.examineRegion.assert_called_once()
    self_mock.selectSingle.assert_not_called()


def test_treated_as_drag_when_x_moves_far():
    self_mock = _make_self_mock()
    canvas = _make_canvas(1000, 600)
    # 500px horizontal drag = 0.5 fraction. The buggy code saw 0.5 < 5 and called
    # selectSingle. The fixed code sees 500px -> examineRegion.
    _call(self_mock, canvas, x0=0.289, y0=0.934, x1=0.746, y1=0.362)
    self_mock.examineRegion.assert_called_once()
    self_mock.selectSingle.assert_not_called()


def test_treated_as_click_on_exact_same_point():
    self_mock = _make_self_mock()
    canvas = _make_canvas(1000, 600)
    _call(self_mock, canvas, x0=0.4, y0=0.4, x1=0.4, y1=0.4, state=1)  # shift-click
    self_mock.selectSingle.assert_called_once()
    self_mock.examineRegion.assert_not_called()


def test_shift_click_drag_to_examineRegion_keeps_state():
    # A shift+drag >= 5px must reach examineRegion with state=shift preserved.
    self_mock = _make_self_mock()
    canvas = _make_canvas(1000, 600)
    _call(self_mock, canvas, x0=0.1, y0=0.1, x1=0.6, y1=0.9, state=1)
    self_mock.examineRegion.assert_called_once()
    assert self_mock.examineRegion.call_args.kwargs["state"] == 1


def test_shift_ctrl_drag_state_not_collapsed_to_no_key():
    # shift+ctrl (1|4 == 5) is the *pick* gesture -> findPeaks inside
    # examineRegion. A port-era "if state != shift and state != ctrl: state =
    # no_key" line folded 5 down to no_key (0), so examineRegion took the
    # search/select branch and pick silently did nothing. The keep-set must
    # include the combined shift+ctrl state. Pin the full 5 survives to the drag
    # path.
    self_mock = _make_self_mock()
    canvas = _make_canvas(1000, 600)
    _call(self_mock, canvas, x0=0.1, y0=0.1, x1=0.6, y1=0.9, state=1 + 4)
    self_mock.examineRegion.assert_called_once()
    assert self_mock.examineRegion.call_args.kwargs["state"] == 1 + 4
    self_mock.selectSingle.assert_not_called()


def test_stray_modifier_bits_are_stripped_to_no_key():
    # On macOS a drag/release's event.state & 255 may carry extra bits a user
    # didn't intend -- CapsLock (0x2), Command (0x10), Option (0x80).  Left
    # alone, those make examineRegion match NO branch, so the box draws but
    # nothing selects (the "box appears but nothing selects" symptom).  The
    # keep-set in selectMulti folds them back to no_key so plain-drag stays
    # robust; pin that.
    for stray in (0x2, 0x10, 0x80, 1 + 0x2, 1 + 4 + 0x2):
        self_mock = _make_self_mock()
        canvas = _make_canvas(1000, 600)
        _call(self_mock, canvas, x0=0.1, y0=0.1, x1=0.6, y1=0.9, state=stray)
        self_mock.examineRegion.assert_called_once()
        assert self_mock.examineRegion.call_args.kwargs["state"] == 0
