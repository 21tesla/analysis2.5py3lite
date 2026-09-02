"""
Unit tests for the macOS 2-finger trackpad scroll -> zoom handler
(WindowFrame.trackpadZoom).

A trackpad scroll arrives as <TouchpadScroll> whose event.delta packs
deltaX (high 16, signed) | deltaY (low 16, signed).  The VERTICAL component
(dy) drives zoom: up-swipe (dy>0) -> zoom in (scale<1), down-swipe (dy<0) ->
zoom out (scale>1), with a continuous scale = exp(-speed * dy).  Modifier
handling mirrors windowsZoom (Ctrl/Shift/Alt scroll the z-plane instead).

These tests drive the unbound method with a mock `self` (only
scrolled_window.zoom / scrollZPlane / touchpad_zoom_speed are needed), so no
tkinter widget and no display are required.
"""

import math
from unittest.mock import MagicMock

from ccpnmr.analysis.frames.WindowFrame import WindowFrame


class _Event:
    def __init__(self, delta, state=0, widget=None):
        self.delta = delta
        self.state = state
        self.widget = widget or MagicMock()


def _self(speed=0.005):
    s = MagicMock()
    s.touchpad_zoom_speed = speed
    return s


def _pack(dx, dy):
    """Build a 32-bit packed delta from signed 16-bit dx, dy (as Tk does)."""
    return ((dx & 0xFFFF) << 16) | (dy & 0xFFFF)


def test_up_swipe_zooms_in_scale_less_than_one():
    self_mock = _self()
    # dy positive -> up -> zoom in (scale < 1)
    event = _Event(_pack(0, 4), state=0)

    WindowFrame.trackpadZoom(self_mock, event)

    self_mock.scrolled_window.zoom.assert_called_once()
    (canvas, scale), _ = self_mock.scrolled_window.zoom.call_args
    assert 0 < scale < 1
    assert math.isclose(scale, math.exp(-0.005 * 4), rel_tol=1e-9)
    self_mock.scrollZPlane.assert_not_called()


def test_down_swipe_zooms_out_scale_greater_than_one():
    self_mock = _self()
    # dy negative -> down -> zoom out (scale > 1)
    event = _Event(_pack(0, -6), state=0)

    WindowFrame.trackpadZoom(self_mock, event)

    (canvas, scale), _ = self_mock.scrolled_window.zoom.call_args
    assert scale > 1
    assert math.isclose(scale, math.exp(-0.005 * -6), rel_tol=1e-9)
    self_mock.scrollZPlane.assert_not_called()


def test_dx_does_not_drive_zoom_only_dy():
    # A large horizontal component with a tiny vertical one must still zoom
    # only by dy (dx is ignored for zoom).
    self_mock = _self(speed=0.01)
    event = _Event(_pack(300, 2), state=0)

    WindowFrame.trackpadZoom(self_mock, event)

    (canvas, scale), _ = self_mock.scrolled_window.zoom.call_args
    # scale depends on dy=2 only
    assert math.isclose(scale, math.exp(-0.01 * 2), rel_tol=1e-9)


def test_zero_dy_is_a_noop():
    self_mock = _self()
    # Pure horizontal scroll (dy == 0): no zoom at all.
    event = _Event(_pack(353, 0), state=0)

    WindowFrame.trackpadZoom(self_mock, event)

    self_mock.scrolled_window.zoom.assert_not_called()
    self_mock.scrollZPlane.assert_not_called()


def test_control_scroll_routes_to_zplane():
    self_mock = _self()
    # Control (0x4) + up swipe -> scrollZPlane(z1, -1), no zoom call.
    event = _Event(_pack(0, 4), state=0x4)

    WindowFrame.trackpadZoom(self_mock, event)

    self_mock.scrolled_window.zoom.assert_not_called()
    self_mock.scrollZPlane.assert_called_once_with(event.widget, "z1", -1)


def test_shift_scroll_routes_to_zplane():
    self_mock = _self()
    # Shift (0x1) + down swipe -> scrollZPlane(z2, +1).
    event = _Event(_pack(0, -4), state=0x1)

    WindowFrame.trackpadZoom(self_mock, event)

    self_mock.scrolled_window.zoom.assert_not_called()
    self_mock.scrollZPlane.assert_called_once_with(event.widget, "z2", 1)


def test_dy_unpacked_as_signed_16bit():
    self_mock = _self(speed=1.0)  # speed=1 so scale == exp(-dy) exactly
    # dy = -1 encoded as 0xFFFF in the low 16 bits must unpack to -1.
    event = _Event(_pack(0, -1), state=0)

    WindowFrame.trackpadZoom(self_mock, event)

    (canvas, scale), _ = self_mock.scrolled_window.zoom.call_args
    # exp(-1 * -1) = exp(1)
    assert math.isclose(scale, math.exp(1.0), rel_tol=1e-9)
