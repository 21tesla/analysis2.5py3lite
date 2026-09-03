# Emulating Linux Input Conventions on macOS Trackpads

> **Status (2026-09-02):** This doc is the *original design exploration* (X11-style
> button emulation, two-finger tap → `Button-2` middle-click, pinch via a
> thumb+index `<B1-Motion>` heuristic). It is **partially superseded** by what was
> actually implemented and verified in `ccpnmr2.5/python/memops/gui/ScrolledWindow.py`
> and `ccpnmr2.5/python/ccpnmr/analysis/frames/WindowFrame.py`. The shipped macOS
> mapping is:
>
> | Gesture | Result | Where |
> | :--- | :--- | :--- |
> | 1-finger tap | point-select peak | `selectSingle` |
> | 1-finger drag | box-select peaks (searchPeaks) | `selectMulti` → `examineRegion` |
> | Ctrl+1-finger drag | box-pick peaks (findPeaks) | `selectMulti` → `examineRegion` |
> | 2-finger **scroll up/down** | zoom (continuous) | `<TouchpadScroll>` → `trackpadZoom` |
> | 2-finger **tap** | open context menu | `macB3*` (release with no commit) |
> | 2-finger **drag** | translate spectrum | `macB3*` (commit to `Button-2` after 5px) |
> | Mouse right **click** | open context menu | `macB3*` (same path as 2-finger tap) |
> | Mouse right **drag** | translate spectrum | `macB3*` (same path as 2-finger drag) |
> | Command / Ctrl / Option + right | open context menu | `macB3*` (on press) |
> | Pinch | **not supported** — no `Magnify`/pinch event reaches Tk | confirmed by probe; 2-finger scroll covers zoom |
>
> Key implementation detail (see `macB3Press/Motion/Release`): a physical right
> click and a trackpad two-finger tap/drag are **byte-identical at the Tk layer**
> (`num=3`, `state=0` on press; no distinguishing field), so the app routes by
> gesture *shape* — a press that moves ≥5px before release is a **drag** (→
> translate); one that stays put is a **tap** (→ menu), decided on release.
> The mapping table in §1 below (e.g. "Control + Tap → `Button-3`") reflects an
> older assumption and no longer drives the code.

---

This document describes the technical architecture, mathematical logic, and edge-case mitigations required to emulate standard Linux/X11 mouse buttons and scrolling behaviors on Apple macOS trackpads (e.g., using Python/Tkinter).

---

## 1. Gesture-to-Button Mapping Matrix

To seamlessly mimic Linux desktop conventions on macOS trackpads, physical trackpad interactions are mapped to virtualized Linux mouse buttons as follows:

| Physical Trackpad Gesture (macOS) | Simulated Linux Event | Description |
| :--- | :--- | :--- |
| **Single-Finger Tap / Click** | `<Button-1>` | Standard Left-Click |
| **Two-Finger Tap / Click** | `<Button-2>` | Middle-Click (Paste, Auto-Scroll) |
| **Control + Tap / Click** | `<Button-3>` | Standard Right-Click |
| **Pinch In (Two Fingers Moving Closer)** | `<Button-4>` | Zoom In |
| **Pinch Out (Two Fingers Moving Apart)** | `<Button-5>` | Zoom Out |
| **Two-Finger Press + Drag** | `<Button-2> Drag` | Middle Drag (Pan / Scroll Lock) |
| **Command + Tap / Click** | `Command-Click` | Standard Log / Utility click (Non-simulated) |

---

## 2. Core Programmatic Logic

### A. Modifier Bitmask Filtering
On macOS, both a **Control-Click** and a physical **Two-Finger Tap** can natively trigger the exact same Right Click event (`<Button-3>` in modern Tkinter engines). To ensure a Two-Finger Tap maps exclusively to `Button-2` (Middle Click) and *never* leaks a `Button-3` event, the state modifier bitmask must be checked:

```python
# Modifier Bitmasks
CONTROL_MASK = 0x4  # State value of 4 represents the physical Control key
COMMAND_MASK = 0x10 # State value of 16 represents the physical Command key

def on_two_finger_press(event):
    # If the Control key modifier is active, this is actually a Control-Click
    if event.state & CONTROL_MASK:
        return "break" # Abort and let the dedicated Control-Click handler run
        
    # Execute pure Two-Finger Click -> Map exclusively to Button-2
    simulate_event("<Button-2>", event)
    return "break" # Prevent event bubbling
```

### B. Two-Finger Vector Pinching
When native trackpad gesture APIs are unavailable or unreliable, pinch-zooming can be emulated by capturing physical thumb-rests + index dragging on a shared `<B1-Motion>` drag stream.

1. **Calculate Distance Vector:**
   On each drag update, compute the Euclidean distance between the anchor point (initial click coordinate, $(x_0, y_0)$) and the moving index finger coordinate $(x_1, y_1)$:
   $$d = \sqrt{(x_1 - x_0)^2 + (y_1 - y_0)^2}$$

2. **Track Change Delta:**
   $$\Delta = d_{\text{current}} - d_{\text{previous}}$$
   * $\Delta < 0$: Fingers are moving closer together (**Pinch In**). Map to `<Button-4>` (Zoom In).
   * $\Delta > 0$: Fingers are moving further apart (**Pinch Out**). Map to `<Button-5>` (Zoom Out).

---

## 3. Special Programmatic Considerations

When writing high-fidelity event emulators, three core engineering challenges must be managed:

### I. Infinite Recursion Loops (Simulation Guards)
When you programmatically generate a virtual mouse event (e.g., `canvas.event_generate("<Button-2>")`), the underlying framework will trigger any bindings listening for that button. If your listener is also the one generating the event, you will trigger an infinite recursion loop, crashing the event loop.

* **Mitigation:** Use localized simulation boolean locks (`_simulating_2` and `_simulating_3`):
  ```python
  def on_two_finger_press(self, event):
      if getattr(self, '_simulating_2', False): 
          return # Abort if we generated this event ourselves
          
      self._simulating_2 = True
      self.canvas.event_generate("<Button-2>", x=event.x, y=event.y)
      self._simulating_2 = False
      return "break"
  ```

### II. The Zoom-Drag Mixing Problem (Gesture State Lock)
Because pinch-zooming and standard clicking/dragging share the exact same trackpad motion stream (`<B1-Motion>`), micro-movement frames and calibration errors will leak standard drag events into your zoom logs.

* **The Problem:**
  1. *Calibration Lag:* On the very first frame of a pinch, no `previous_distance` exists yet, causing the code to default to logging a standard drag coordinate.
  2. *Jitter Thresholds:* High-frequency trackpad polling means some updates have a tiny movement delta ($\le 5\text{px}$). Standard filtering skips these, causing the event to "fall through" and trigger standard drag actions.

* **The Mitigation (The `is_pinching` State Lock):**
  Initialize a boolean flag `is_pinching = False` on click. Once the finger distance delta crosses the initial calibration threshold ($> 5\text{px}$), lock the state to `True`. 
  
  As long as the state is locked, **all standard drag operations are entirely bypassed**, regardless of whether subsequent frame deltas drop below the threshold:

  ```python
  def on_drag_pipeline(self, event):
      # Skip standard drag paths if we are locked into a zoom operation
      if self.is_pinching or abs(delta) > 5:
          self.is_pinching = True # Lock the state
          
          if abs(delta) > 5:
              # Emit Button-4/Button-5 event
              ...
          return # ALWAYS exit here to prevent standard drag leaking
          
      # Run standard dragging only if the lock is False
      log_standard_drag(event)
  ```

### III. Event Swallowing (`return "break"`)
In GUI frameworks like Tkinter, events bubble up through class and window hierarchies. If you translate a Control-Click to a `<Button-3>` but do not stop the original event, both the Control-Click and standard click handlers may fire.
* **Mitigation:** Always return `"break"` at the end of custom mapped handlers to indicate that the event is fully consumed and should not propagate to other bound listeners.
