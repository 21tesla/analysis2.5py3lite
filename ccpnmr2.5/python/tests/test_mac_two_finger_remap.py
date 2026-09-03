"""
Unit tests for the macOS Button-3 (right button / two-finger) tap-vs-drag
state machine (ScrolledWindow.macB3Press / macB3Motion / macB3Release).

On macOS a physical right-click and a trackpad two-finger press are
byte-identical at the Tk layer (probe btn3_probe.py, 2026-09-02, Tk 9.0.4:
press num=3 state=0; release num=3 state=1024 = button mask; no keycode/char/
delta distinguish them, for either a bare tap or a drag). So the app routes by
GESTURE SHAPE, not source, via a deferred state machine:
  - no-modifier Button-3 press  -> remember origin, wait for motion/release
  - first B3-Motion past a 5px threshold -> commit to a drag: fire <ButtonPress-2>
    and forward motions to <B2-Motion> (the app's Button-2 translate pipeline)
  - Button-3 release with no commit -> a TAP -> open the action menu
  - Command/Control/Option press -> open the action menu immediately

These tests drive the (unbound by default) methods with a mock `self`, so no
tkinter widget or display is required.
"""

from unittest.mock import MagicMock

from memops.gui.ScrolledWindow import ScrolledWindow

# Button-2 state mask used for a generated <B2-Motion>.
B2MASK = 1 << (7 + 2)
# Matches ScrolledWindow._mac_b3_drag_threshold.
THRESHOLD = 5


class _Event:
    def __init__(self, state, x=55, y=60, widget=None):
        self.state = state
        self.x = x
        self.y = y
        self.widget = widget or MagicMock()


def _self():
    m = MagicMock()
    # Give the state fields honest starting values (a fresh ScrolledWindow).
    m._mac_b3_pressed = False
    m._mac_b3_x0 = 0
    m._mac_b3_y0 = 0
    m._mac_two_finger = False
    # The threshold is a CLASS attribute on ScrolledWindow, but the test passes
    # a bare MagicMock as `self`, so seed the real value here.
    m._mac_b3_drag_threshold = ScrolledWindow._mac_b3_drag_threshold
    return m


def test_no_modifier_press_remembers_origin_and_does_nothing_else():
    self_mock = _self()
    event = _Event(state=0, x=100, y=120)

    ret = ScrolledWindow.macB3Press(self_mock, event)

    # No menu, no <ButtonPress-2>: the tap-vs-drag outcome is deferred.
    self_mock.menu.popupMenu.assert_not_called()
    event.widget.event_generate.assert_not_called()
    assert self_mock._mac_b3_pressed is True
    assert self_mock._mac_b3_x0 == 100
    assert self_mock._mac_b3_y0 == 120
    assert self_mock._mac_two_finger is False
    assert ret == "break"


def test_motion_within_threshold_is_ignored():
    self_mock = _self()
    ScrolledWindow.macB3Press(self_mock, _Event(state=0, x=100, y=120))
    # 3px in x, 0px in y: below the 5px threshold -> still a tap candidate.
    event = _Event(state=B2MASK | 1 << (7 + 3), x=103, y=120)

    ret = ScrolledWindow.macB3Motion(self_mock, event)

    event.widget.event_generate.assert_not_called()
    assert self_mock._mac_two_finger is False
    assert ret == "break"


def test_motion_past_threshold_commits_to_b2_translate():
    self_mock = _self()
    ScrolledWindow.macB3Press(self_mock, _Event(state=0, x=100, y=120))
    # 10px in x: crosses the threshold.
    event = _Event(state=B2MASK | 1 << (7 + 3), x=110, y=120)

    ret = ScrolledWindow.macB3Motion(self_mock, event)

    gens = [g.args[0] for g in event.widget.event_generate.call_args_list]
    # ButtonPress-2 is fired first (at the press origin), then the B2-Motion.
    assert gens == ["<ButtonPress-2>", "<B2-Motion>"]
    press_gen = event.widget.event_generate.call_args_list[0]
    assert press_gen.kwargs.get("button") == 2
    assert press_gen.kwargs.get("x") == 100  # origin, not the current x
    assert press_gen.kwargs.get("y") == 120
    move_gen = event.widget.event_generate.call_args_list[1]
    assert (move_gen.kwargs.get("state") & B2MASK) == B2MASK  # carries Button-2 mask
    assert move_gen.kwargs.get("x") == 110
    assert move_gen.kwargs.get("y") == 120
    assert self_mock._mac_two_finger is True
    assert ret == "break"


def test_subsequent_motion_after_commit_forwards_b2_only():
    self_mock = _self()
    ScrolledWindow.macB3Press(self_mock, _Event(state=0, x=100, y=120))
    ScrolledWindow.macB3Motion(self_mock, _Event(state=B2MASK | 1 << (7 + 3), x=120, y=120))  # commits
    event = _Event(state=B2MASK | 1 << (7 + 3), x=130, y=125)

    ret = ScrolledWindow.macB3Motion(self_mock, event)

    # Already committed: only the B2-Motion fires, no second ButtonPress-2.
    gens = [g.args[0] for g in event.widget.event_generate.call_args_list]
    assert gens == ["<B2-Motion>"]
    assert (event.widget.event_generate.call_args.kwargs.get("state") & B2MASK) == B2MASK
    assert ret == "break"


def test_motion_without_press_is_noop():
    self_mock = _self()
    event = _Event(state=B2MASK | 1 << (7 + 3), x=10, y=10)

    ret = ScrolledWindow.macB3Motion(self_mock, event)

    event.widget.event_generate.assert_not_called()
    assert ret == "break"


def test_release_after_commit_closes_b2_sequence():
    self_mock = _self()
    ScrolledWindow.macB3Press(self_mock, _Event(state=0, x=100, y=120))
    ScrolledWindow.macB3Motion(self_mock, _Event(state=B2MASK | 1 << (7 + 3), x=120, y=120))  # commits
    event = _Event(state=0, x=140, y=130)

    ret = ScrolledWindow.macB3Release(self_mock, event)

    gens = [g.args[0] for g in event.widget.event_generate.call_args_list]
    assert gens == ["<ButtonRelease-2>"]
    assert event.widget.event_generate.call_args.kwargs.get("button") == 2
    self_mock.menu.popupMenu.assert_not_called()
    assert self_mock._mac_two_finger is False
    assert ret == "break"


def test_release_with_no_commit_is_a_tap_that_opens_menu():
    self_mock = _self()
    ScrolledWindow.macB3Press(self_mock, _Event(state=0, x=100, y=120))
    event = _Event(state=0, x=100, y=120)

    ret = ScrolledWindow.macB3Release(self_mock, event)

    # A stationary tap opens the menu.
    self_mock.menu.popupMenu.assert_called_once_with(event)
    event.widget.event_generate.assert_not_called()
    assert self_mock._mac_b3_pressed is False
    assert self_mock._mac_two_finger is False
    assert ret == "break"


def test_release_after_small_motion_is_still_a_tap():
    self_mock = _self()
    ScrolledWindow.macB3Press(self_mock, _Event(state=0, x=100, y=120))
    # 2px jiggle, below threshold: no commit.
    ScrolledWindow.macB3Motion(self_mock, _Event(state=B2MASK | 1 << (7 + 3), x=102, y=120))
    event = _Event(state=0, x=102, y=120)

    ret = ScrolledWindow.macB3Release(self_mock, event)

    self_mock.menu.popupMenu.assert_called_once_with(event)
    event.widget.event_generate.assert_not_called()
    assert self_mock._mac_two_finger is False
    assert ret == "break"


def test_release_with_no_state_is_noop():
    self_mock = _self()
    event = _Event(state=0)

    ret = ScrolledWindow.macB3Release(self_mock, event)

    event.widget.event_generate.assert_not_called()
    self_mock.menu.popupMenu.assert_not_called()
    assert ret == "break"


def test_control_click_opens_menu_immediately():
    self_mock = _self()
    event = _Event(state=0x4, x=12, y=34)  # Control bit (& 255 -> 4)

    ret = ScrolledWindow.macB3Press(self_mock, event)

    self_mock.menu.popupMenu.assert_called_once_with(event)
    event.widget.event_generate.assert_not_called()
    assert self_mock._mac_b3_pressed is False
    assert self_mock._mac_two_finger is False
    assert ret == "break"


def test_control_click_with_button_mask_bits():
    self_mock = _self()
    # Control (0x4) plus the Button-3 mask bits Tk adds: & 255 keeps the Control.
    event = _Event(state=0x4 | (1 << (7 + 3)))

    ScrolledWindow.macB3Press(self_mock, event)

    self_mock.menu.popupMenu.assert_called_once_with(event)
    event.widget.event_generate.assert_not_called()


def test_command_click_opens_menu_immediately():
    self_mock = _self()
    event = _Event(state=0x10)  # Command

    ScrolledWindow.macB3Press(self_mock, event)

    self_mock.menu.popupMenu.assert_called_once_with(event)
    event.widget.event_generate.assert_not_called()
    assert self_mock._mac_b3_pressed is False
    assert self_mock._mac_two_finger is False


def test_option_click_opens_menu_immediately():
    self_mock = _self()
    event = _Event(state=0x80)  # Option

    ScrolledWindow.macB3Press(self_mock, event)

    self_mock.menu.popupMenu.assert_called_once_with(event)
    event.widget.event_generate.assert_not_called()
    assert self_mock._mac_b3_pressed is False


def test_shift_only_press_is_not_menu_but_is_tracked():
    self_mock = _self()
    # Shift is not a menu modifier: it should arm the tap-vs-drag path, NOT open
    # the menu, so a shift+right tap/drag still resolves by shape.
    event = _Event(state=0x1, x=50, y=60)

    ret = ScrolledWindow.macB3Press(self_mock, event)

    self_mock.menu.popupMenu.assert_not_called()
    assert self_mock._mac_b3_pressed is True
    assert self_mock._mac_b3_x0 == 50
    assert ret == "break"
