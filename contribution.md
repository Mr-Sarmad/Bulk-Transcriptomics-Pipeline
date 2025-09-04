# Contributing Guidelines

Thank you for your interest in contributing to the **RNA-seq Analysis Pipeline**  
This project aims to provide a simple, reproducible, and computationally friendly RNA-seq workflow that can run even on minimal hardware setups.

---

## How to Contribute

1. **Fork the Repository**
   - Click the **Fork** button at the top-right of this page.
   - Clone your fork to your local machine:
     ```bash
     git clone https://github.com/Mr-Sarmad/rnaseq-pipeline.git
     cd rnaseq-pipeline
     ```

2. **Set Up Your Environment**
   - Create a conda environment using the provided file:
     ```bash
     conda env create -f environment.yaml
     conda activate rnaseq
     ```
   - Or install requirements directly:
     ```bash
     pip install -r requirements.txt
     ```

3. **Make Your Changes**
   - Work on a new branch:
     ```bash
     git checkout -b feature/my-feature
     ```
   - Make your improvements (e.g., new analysis step, bug fix, documentation).

4. **Test the Pipeline**
   - Ensure the Snakemake workflow runs without errors:
     ```bash
     snakemake --cores 4
     ```

5. **Submit a Pull Request**
   - Push your branch:
     ```bash
     git push origin feature/my-feature
     ```
   - Open a Pull Request (PR) on GitHub and describe your changes.

---

## Reporting Issues

- If you find a bug, please open an issue in the **Issues** tab.  
- Provide:
  - A clear description of the problem  
  - Steps to reproduce it  
  - Error messages or logs (if available)

---

## Code Style

- Keep scripts clear and modular (R scripts go in `scripts/` if added later).  
- Use informative commit messages (e.g., `fix: corrected tximport path handling`).  
- Follow existing file structure:
