import os
import sys
import subprocess
import webbrowser
import time

def main():
    print("🚀 Starting Aarvi Chem AI server...")

    # current folder
    base_dir = os.path.dirname(os.path.abspath(__file__))

    # run uvicorn
    cmd = [
        sys.executable,
        "-m",
        "uvicorn",
        "api:app",
        "--host", "127.0.0.1",
        "--port", "8000"
    ]

    subprocess.Popen(cmd, cwd=base_dir)

    # wait a bit
    time.sleep(3)

    # open browser
    webbrowser.open("http://127.0.0.1:8000")

if __name__ == "__main__":
    main()
