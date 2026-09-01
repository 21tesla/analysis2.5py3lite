"""Unit tests for WindowFrame mouse wheel zoom and Z plane scrolling bindings and functions.
"""
import pytest
from unittest.mock import MagicMock, patch

from ccpnmr.analysis.frames.WindowFrame import WindowFrame


class MockEvent:
    def __init__(self, delta, state, widget=None):
        self.delta = delta
        self.state = state
        self.widget = widget or MagicMock()


def test_windows_zoom_no_modifiers():
    # Mock self/WindowFrame
    self_mock = MagicMock(spec=WindowFrame)
    self_mock.scrolled_window = MagicMock()
    self_mock.scrollZPlane = MagicMock()

    # Create dummy canvas widget
    canvas = MagicMock()

    # 1. Test Zoom In (delta > 0)
    event_in = MockEvent(delta=120, state=0, widget=canvas)
    WindowFrame.windowsZoom(self_mock, event_in)

    # Assert zoomed with scale 0.8
    self_mock.scrolled_window.zoom.assert_called_once_with(canvas, 0.8)
    self_mock.scrollZPlane.assert_not_called()

    # Reset mocks
    self_mock.scrolled_window.zoom.reset_mock()
    self_mock.scrollZPlane.reset_mock()

    # 2. Test Zoom Out (delta < 0)
    event_out = MockEvent(delta=-120, state=0, widget=canvas)
    WindowFrame.windowsZoom(self_mock, event_out)

    # Assert zoomed with scale 1.2
    self_mock.scrolled_window.zoom.assert_called_once_with(canvas, 1.2)
    self_mock.scrollZPlane.assert_not_called()


def test_windows_zoom_control_modifier():
    self_mock = MagicMock(spec=WindowFrame)
    self_mock.scrolled_window = MagicMock()
    self_mock.scrollZPlane = MagicMock()

    canvas = MagicMock()

    # Control modifier (state & 4) - Zoom In (delta > 0)
    event_ctrl_in = MockEvent(delta=120, state=4, widget=canvas)
    WindowFrame.windowsZoom(self_mock, event_ctrl_in)

    # Assert scrolled z1 plane with step -1
    self_mock.scrollZPlane.assert_called_once_with(canvas, "z1", -1)
    self_mock.scrolled_window.zoom.assert_not_called()

    # Reset mocks
    self_mock.scrollZPlane.reset_mock()

    # Control modifier (state & 4) - Zoom Out (delta < 0)
    event_ctrl_out = MockEvent(delta=-120, state=4, widget=canvas)
    WindowFrame.windowsZoom(self_mock, event_ctrl_out)

    # Assert scrolled z1 plane with step 1
    self_mock.scrollZPlane.assert_called_once_with(canvas, "z1", 1)
    self_mock.scrolled_window.zoom.assert_not_called()


def test_windows_zoom_shift_modifier():
    self_mock = MagicMock(spec=WindowFrame)
    self_mock.scrolled_window = MagicMock()
    self_mock.scrollZPlane = MagicMock()

    canvas = MagicMock()

    # Shift modifier (state & 1) - Zoom In (delta > 0)
    event_shift_in = MockEvent(delta=120, state=1, widget=canvas)
    WindowFrame.windowsZoom(self_mock, event_shift_in)

    # Assert scrolled z2 plane with step -1
    self_mock.scrollZPlane.assert_called_once_with(canvas, "z2", -1)
    self_mock.scrolled_window.zoom.assert_not_called()

    # Reset mocks
    self_mock.scrollZPlane.reset_mock()

    # Shift modifier (state & 1) - Zoom Out (delta < 0)
    event_shift_out = MockEvent(delta=-120, state=1, widget=canvas)
    WindowFrame.windowsZoom(self_mock, event_shift_out)

    # Assert scrolled z2 plane with step 1
    self_mock.scrollZPlane.assert_called_once_with(canvas, "z2", 1)
    self_mock.scrolled_window.zoom.assert_not_called()


def test_windows_zoom_alt_modifier():
    self_mock = MagicMock(spec=WindowFrame)
    self_mock.scrolled_window = MagicMock()
    self_mock.scrollZPlane = MagicMock()

    canvas = MagicMock()

    # Alt modifier (state & 8) - Zoom In (delta > 0)
    event_alt_in = MockEvent(delta=120, state=8, widget=canvas)
    WindowFrame.windowsZoom(self_mock, event_alt_in)

    # Assert scrolled z1 plane with step -1 (Alt scrolls z1)
    self_mock.scrollZPlane.assert_called_once_with(canvas, "z1", -1)
    self_mock.scrolled_window.zoom.assert_not_called()

    # Reset mocks
    self_mock.scrollZPlane.reset_mock()

    # Alt modifier (state & 8) - Zoom Out (delta < 0)
    event_alt_out = MockEvent(delta=-120, state=8, widget=canvas)
    WindowFrame.windowsZoom(self_mock, event_alt_out)

    # Assert scrolled z1 plane with step 1
    self_mock.scrollZPlane.assert_called_once_with(canvas, "z1", 1)
    self_mock.scrolled_window.zoom.assert_not_called()


def test_windows_zoom_zero_delta():
    self_mock = MagicMock(spec=WindowFrame)
    self_mock.scrolled_window = MagicMock()
    self_mock.scrollZPlane = MagicMock()

    canvas = MagicMock()

    event_zero = MockEvent(delta=0, state=0, widget=canvas)
    WindowFrame.windowsZoom(self_mock, event_zero)

    # Assert neither scroll nor zoom is called
    self_mock.scrollZPlane.assert_not_called()
    self_mock.scrolled_window.zoom.assert_not_called()
