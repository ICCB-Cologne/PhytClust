# Phytclust Web Frontend
This Flask app is a prototype of a PhytClust web frontend.
Run with `python3 app.py`. Paste a valid Newick string into
the text entry field to see the tree visualized and a
textual representation of the best clusters that PhytClust
finds.


# Installation

Install local version of PhytClust with `pip install -e .` and install dependencies using `pip install "uvicorn[standard]" fastapi`.


# Running on local

Start local server with `uvicorn phytclust.gui.api:app --reload --host 127.0.0.1 --port 8000`. Then open `http://127.0.0.1:8000` in a browser.


