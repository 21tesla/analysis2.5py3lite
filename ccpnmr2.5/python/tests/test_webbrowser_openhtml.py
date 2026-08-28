"""Project > Summary -> WebBrowser.openHtml regression.

Symptom (the user's error.txt): clicking Project > Summary raised
``TypeError: a bytes-like object is required, not 'str'`` in
``WebBrowser.openHtml`` — ``tempfile.NamedTemporaryFile()`` defaults to
binary mode, so the str write dies on py3.  The fix must also keep the
file at its path AFTER openHtml returns (the browser follows the file://
URL asynchronously, so the default delete-on-create unlink would leave it
dangling).
"""
import os


class TestWebBrowserOpenHtml:
    def test_str_writes_to_file_still_readable_after_return(self):
        from memops.gui import WebBrowser

        opened = []
        wb = WebBrowser.WebBrowser(parent=None, name="firefox")
        wb.open = lambda url: opened.append(url)

        html = "<html><body><table><td>hello</td></table></body></html>"
        wb.openHtml(html)  # raised TypeError: a bytes-like object ... pre-fix

        assert len(opened) == 1
        url = opened[0]
        assert url.startswith("file://")
        path = url[len("file://"):]
        try:
            assert os.path.exists(path)  # delete=False: browser reads it after close
            assert path.endswith(".html")  # browsers sniff the type from the extension
            with open(path) as f:
                assert f.read() == html
        finally:
            if os.path.exists(path):
                os.unlink(path)
