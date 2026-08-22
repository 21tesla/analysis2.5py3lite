import tkinter as Tkinter

import time
import sys

def LOG(msg):
    with open('/tmp/ccpn_modal_debug.log', 'a') as f:
        f.write(f"[{time.time():.3f}] {msg}\n")
    print(f"DEBUG: {msg}")
    sys.stdout.flush()

from memops.gui.Base import Base
from memops.gui.ButtonList import ButtonList
from memops.gui.Entry import Entry
from memops.gui.Label import Label
from memops.gui.MessageReporter import showWarning


class QueryDialogBox(Tkinter.Toplevel, Base):
    def __init__(
        self,
        parent,
        title="Dialog Box",
        prompt="Prompt",
        initialvalue=None,
        minvalue=None,
        maxvalue=None,
        position=(50, 50),
        returnType=str,
        hide=0,
    ):
        LOG(f"QueryDialogBox.__init__ started. title={title}, parent={parent}")
        if parent is None:
            self.parent = Tkinter._default_root
        else:
            self.parent = parent

        LOG("Calling Tkinter.Toplevel.__init__")
        Tkinter.Toplevel.__init__(self, parent)
        LOG("Calling self.withdraw()")
        self.withdraw() # Start hidden to prevent macOS rendering races

        self.response = None
        self.prompt = prompt
        self.initVal = initialvalue
        self.minVal = minvalue
        self.maxVal = maxvalue
        self.hide = hide
        self.returnType = returnType

        LOG("Calling self.transient()")
        self.transient(parent)
        LOG("Calling self.title()")
        self.title(title or "")

        LOG("Creating mainFrame")
        mainFrame = Tkinter.Frame(self)
        LOG("Packing mainFrame")
        mainFrame.pack(padx=2, pady=2)

        LOG("Calling drawButtons")
        self.drawButtons("OK", "Cancel")
        LOG("Calling body")
        self.entry = self.body(mainFrame)

        LOG("Calling protocol and bind")
        self.protocol("WM_WINDOW_DELETE", self.cancel)
        self.bind("<Return>", self.ok)
        self.bind("<Escape>", self.cancel)

        if parent is not None:
            LOG("Setting geometry")
            dx = parent.winfo_rootx() + position[0]
            dy = parent.winfo_rooty() + position[1]
            self.geometry("+%d+%d" % (dx, dy))

        LOG("Calling focus_set on entry")
        self.entry.focus_set()

        LOG("Calling deiconify")
        self.deiconify() # Display window
        LOG("Calling update_idletasks")
        self.update_idletasks() # Flush rendering
        LOG("Calling wait_visibility")
        self.wait_visibility() # Wait for OS confirmation
        LOG("Calling grab_set")
        self.grab_set()  # Make modal safely

        LOG("Setting up wait_variable")
        self.var = Tkinter.IntVar()
        LOG("Calling wait_variable")
        self.wait_variable(self.var)
        LOG("wait_variable finished!")

    def body(self, frame):
        LOG("body: Creating label")
        label = Label(frame, text=self.prompt)

        LOG("body: Creating entry")
        if self.hide:
            entry = Entry(frame, text=self.initVal, bg="white", show="*")
        else:
            entry = Entry(frame, text=self.initVal, bg="white")

        LOG("body: grid label")
        label.grid(row=0, column=0, sticky=Tkinter.W)
        LOG("body: grid entry")
        entry.grid(row=1, column=0, sticky=Tkinter.W)

        LOG("body: calling focus_set")
        entry.focus_set()  # for the keypresses
        LOG("body: returning entry")
        return entry

    def cancel(self, *event):
        LOG("cancel() called")
        if self.parent:
            self.parent.focus_set()
        if hasattr(self, 'var'):
            self.var.set(1)
        Tkinter.Toplevel.destroy(self)
        LOG("cancel() finished")

    def ok(self, *event):
        LOG("ok() called")
        self.response = self.getResponse()
        if self.response is not None:
            LOG(f"Response is {self.response}")
            self.update_idletasks()
            if hasattr(self, 'var'):
                self.var.set(1)
            self.destroy()
            LOG("ok() destroyed self")
        else:
            LOG("Response is None")
            self.parent.focus_set()

    def drawButtons(self, okText="OK", cancelText="Cancel"):
        LOG("drawButtons: starting")
        texts = [okText, cancelText]
        commands = [self.ok, self.cancel]

        LOG("drawButtons: creating ButtonList")
        buttonList = ButtonList(self, texts=texts, commands=commands)
        LOG("drawButtons: configuring default button")
        buttonList.buttons[0].config(default=Tkinter.ACTIVE)
        LOG("drawButtons: packing buttonList")
        buttonList.pack()
        LOG("drawButtons: finished")

    def getResponse(self):

        response = self.entry.get()
        if self.returnType == int:
            try:
                response = int(response)
            except:
                raise Exception("Value %s not an integer" % response)

            response = self.checkValue(response)

        elif self.returnType == float:
            try:
                response = float(response)
            except:
                raise Exception("Value %s not a floating point number" % response)

            response = self.checkValue(response)

        elif self.returnType == str:
            response = response.strip()

        else:
            raise Exception("Unknown return type for query box")

        return response

    def checkValue(self, value):

        if (self.minVal is not None) and (value < self.minVal):
            showWarning("Warning", "Value is less than minimum value (%s)" % self.minVal)
            value = None

        elif (self.maxVal is not None) and (value > self.maxVal):
            showWarning("Warning", "Value is greater than maximum value (%s)" % self.maxVal)
            value = None

        return value


def askInteger(title, prompt, parent=None, **kw):
    kw["returnType"] = int
    dialog = QueryDialogBox(parent, title, prompt, **kw)
    return dialog.response


def askFloat(title, prompt, parent=None, **kw):
    kw["returnType"] = float
    dialog = QueryDialogBox(parent, title, prompt, **kw)
    return dialog.response


def askString(title, prompt, parent=None, **kw):
    kw["returnType"] = str
    dialog = QueryDialogBox(parent, title, prompt, **kw)
    return dialog.response


def askPassword(title, prompt, parent=None, **kw):
    kw["returnType"] = str
    dialog = QueryDialogBox(parent, title, prompt, hide=1, **kw)
    return dialog.response
