"""
Unit tests for the macOS two-finger -> middle (Button-2) remap
(ScrolledWindow.macTwoFingerToMiddle / macTwoFingerMotion / macTwoFingerRelease).

On macOS a trackpad two-finger tap and a physical right-click both arrive
as <ButtonPress-3>, and a two-finger tap-drag arrives as
<ButtonPress-3> + <B3-Motion>* + <ButtonRelease-3>.  The remap distinguishes
by modifier:
  - Control / Command / Option held  -> open the action menu (right-click)
  - no such modifier -> remap the whole Button-3 sequence to a middle
    (Button-2) sequence so the app's normal Button-2 handlers run: a tap
    becomes a neutral middle click (press -> release, no translate) and a
    tap-drag becomes a middle drag that translates the spectrum.

These tests drive the unbound methods with a mock `self` that only provides
`menu` and `_mac_two_finger`, so no tkinter widget (and no display) is
required.  (The full press->motion->release flow is exercised against a real
Tk canvas in manual verification.)
"""

from unittest.mock import MagicMock

from memops.gui.ScrolledWindow import ScrolledWindow

# Button-2 state mask used for a generated <B2-Motion>
B2MASK = 1 << (7 + 2)


class _Event:
    def __init__(self, state, x=55, y=60, widget=None):
        self.state = state
        self.x = x
        self.y = y
        self.widget = widget or MagicMock()


def _self():
    return MagicMock()


def test_two_finger_press_generates_middle_press():
    self_mock = _self()
    event = _Event(state=0)  # no modifiers: bare two-finger press

    ret = ScrolledWindow.macTwoFingerToMiddle(self_mock, event)

    # Action menu must NOT open on a two-finger press.
    self_mock.menu.popupMenu.assert_not_called()
    # Only a middle Press-2 is generated here; the release is produced by
    # macTwoFingerRelease when the user lets go.
    gens = event.widget.event_generate.call_args_list
    assert [g.args[0] for g in gens] == ["<ButtonPress-2>"]
    for g in gens:
        assert g.kwargs.get("x") == 55
        assert g.kwargs.get("y") == 60
        assert g.kwargs.get("button") == 2
    # Remap is armed for the follow-on motion/release handlers.
    assert self_mock._mac_two_finger is True
    assert ret == "break"


def test_two_finger_release_closes_middle_and_disarms():
    self_mock = _self()
    self_mock._mac_two_finger = True
    event = _Event(state=0, x=70, y=80)

    ret = ScrolledWindow.macTwoFingerRelease(self_mock, event)

    assert [g.args[0] for g in event.widget.event_generate.call_args_list] == ["<ButtonRelease-2>"]
    assert event.widget.event_generate.call_args.kwargs["button"] == 2
    assert self_mock._mac_two_finger is False
    assert ret == "break"


def test_two_finger_release_is_noop_when_not_remapping():
    self_mock = _self()
    self_mock._mac_two_finger = False
    event = _Event(state=0)

    ret = ScrolledWindow.macTwoFingerRelease(self_mock, event)

    event.widget.event_generate.assert_not_called()
    assert self_mock._mac_two_finger is False
    assert ret == "break"


def test_two_finger_motion_generates_middle_motion_with_mask():
    self_mock = _self()
    self_mock._mac_two_finger = True
    event = _Event(state=1 << (7 + 3), x=60, y=72)  # B3-Motion carries its B3 mask

    ret = ScrolledWindow.macTwoFingerMotion(self_mock, event)

    gens = event.widget.event_generate.call_args_list
    assert [g.args[0] for g in gens] == ["<B2-Motion>"]
    kwargs = gens[0].kwargs
    # A generated <B2-Motion> needs the button-2 mask in -state to be
    # recognised as a button-2 motion (state 0 demotes to plain <Motion>).
    assert (kwargs.get("state") & B2MASK) == B2MASK
    assert kwargs.get("x") == 60
    assert kwargs.get("y") == 72
    assert self_mock._mac_two_finger is True
    assert ret == "break"


def test_two_finger_motion_is_noop_when_not_remapping():
    self_mock = _self()
    self_mock._mac_two_finger = False
    event = _Event(state=1 << (7 + 3))

    ret = ScrolledWindow.macTwoFingerMotion(self_mock, event)

    event.widget.event_generate.assert_not_called()
    assert ret == "break"


def test_control_click_opens_menu_not_middle():
    self_mock = _self()
    # Real Tk events carry the button mask high bits too; the handler must
    # mask with & 255 so a Control press is still recognised as Control.
    event = _Event(state=0x4, x=12, y=34)  # Control

    ret = ScrolledWindow.macTwoFingerToMiddle(self_mock, event)

    # Menu requested with the original Button-3 event.
    self_mock.menu.popupMenu.assert_called_once_with(event)
    # No middle (Button-2) event is generated for a Control click, and no
    # remap is armed (clears any stale state).
    event.widget.event_generate.assert_not_called()
    assert self_mock._mac_two_finger is False
    assert ret == "break"


def test_control_click_with_button_mask_bits():
    self_mock = _self()
    # Control (0x4) plus the Button-3 mask bits Tk adds on macOS.  After
    # masking with & 255 the Control bit must remain set -> menu.
    event = _Event(state=0x4 | 0x200 | 0x800)

    ScrolledWindow.macTwoFingerToMiddle(self_mock, event)

    self_mock.menu.popupMenu.assert_called_once_with(event)
    event.widget.event_generate.assert_not_called()


def test_command_click_opens_menu_not_middle():
    self_mock = _self()
    event = _Event(state=0x10)  # Command

    ScrolledWindow.macTwoFingerToMiddle(self_mock, event)

    self_mock.menu.popupMenu.assert_called_once_with(event)
    event.widget.event_generate.assert_not_called()
    assert self_mock._mac_two_finger is False


def test_option_click_opens_menu_not_middle():
    self_mock = _self()
    event = _Event(state=0x80)  # Option

    ScrolledWindow.macTwoFingerToMiddle(self_mock, event)

    self_mock.menu.popupMenu.assert_called_once_with(event)
    event.widget.event_generate.assert_not_called()
    assert self_mock._mac_two_finger is False


def test_shift_only_press_remapped_not_menu():
    self_mock = _self()
    # Shift alone is not a menu modifier -> treated as a (shift) two-finger
    # press, i.e. remapped to middle, not the menu.
    event = _Event(state=0x1)

    ret = ScrolledWindow.macTwoFingerToMiddle(self_mock, event)

    self_mock.menu.popupMenu.assert_not_called()
    assert [g.args[0] for g in event.widget.event_generate.call_args_list] == ["<ButtonPress-2>"]
    assert self_mock._mac_two_finger is True
    assert ret == "break"
