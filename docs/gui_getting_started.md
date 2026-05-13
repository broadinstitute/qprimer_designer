# Getting Started with the qPrimer Designer GUI

This guide walks you through downloading the code from GitHub and launching the graphical interface on **Mac** or **Windows**.

---

## Prerequisites

You need two things installed before starting:

1. **Git** -- to download the code
2. **Conda** (via Miniforge) -- to install Python and dependencies

### Install Git

**Mac:** Open **Terminal** (search "Terminal" in Spotlight). Type `git --version`. If not installed, macOS will prompt you to install the Command Line Tools -- follow the prompts.

**Windows:** Download and install [Git for Windows](https://gitforwindows.org/). Use all default settings. After installation, open **Git Bash** from the Start menu.

### Install Conda (Miniforge)

Miniforge is a lightweight installer for conda. Download it from [github.com/conda-forge/miniforge](https://github.com/conda-forge/miniforge#miniforge3).

**Mac:** Choose the macOS **arm64 (Apple Silicon)** or **x86_64 (Intel)** version depending on your Mac. Open Terminal and run: `bash Miniforge3-MacOSX-arm64.sh` (or the Intel filename). Follow the prompts, then restart Terminal.

**Windows:** Download the **Windows x86_64** installer from the same page. Double-click the `.exe` and follow the prompts. After installation, open **Miniforge Prompt** from the Start menu.

Verify conda is installed by typing:
```bash
conda --version
```

---

## Step 1: Download the code

Open your terminal (**Terminal** on Mac, **Miniforge Prompt** on Windows) and run:

```bash
git clone https://github.com/broadinstitute/qprimer_designer.git
cd qprimer_designer
```

This downloads the entire project into a folder called `qprimer_designer`.

---

## Step 2: Create the environment and install

```bash
conda env create -f environment.yml
conda activate qprimer-designer
pip install -e ".[gui]"
```

This installs Python, all scientific dependencies, and the GUI libraries. It may take several minutes the first time.

---

## Step 3: Launch the GUI

```bash
streamlit run gui/app.py
```

After a few seconds, you should see output like:

```
  You can now view your Streamlit app in your browser.

  Local URL: http://localhost:8501
```

Your web browser should open automatically. If it doesn't, open your browser and go to **http://localhost:8501**.

---

## Stopping the GUI

Go back to the terminal and press **Ctrl+C** to stop the server.

---

## Returning after the first time

After initial setup, you only need two commands to relaunch:

```bash
conda activate qprimer-designer
streamlit run gui/app.py
```

(Run these from inside the `qprimer_designer` folder.)

---

## Troubleshooting

**`conda: command not found`** -- Restart your terminal after installing Miniforge. On Windows, make sure you're using **Miniforge Prompt**, not regular Command Prompt.

**`git: command not found`** -- On Mac, run `xcode-select --install`. On Windows, make sure you installed Git for Windows and are using Git Bash or Miniforge Prompt.

**`streamlit: command not found`** -- Make sure you ran `pip install -e ".[gui]"` and that the conda environment is active (`conda activate qprimer-designer`).

**Browser doesn't open** -- Manually navigate to `http://localhost:8501` in your browser.

**Environment creation is slow** -- This is normal for the first install. Subsequent activations will be fast.
