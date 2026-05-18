#!/usr/bin/env python
"""
setup_and_run.py
----------------
Installation helper and setup script for AI-Powered Plant Gene Analyzer.

This script:
1. Checks Python version
2. Creates virtual environment (optional)
3. Installs dependencies from requirements.txt
4. Runs health checks
5. Launches the Streamlit app

Run with: python setup_and_run.py
"""

import subprocess
import sys
import os
from pathlib import Path


def print_header(text: str):
    """Print a formatted header."""
    print("\n" + "=" * 70)
    print(f"🧬 {text}")
    print("=" * 70)


def check_python_version():
    """Verify Python version >= 3.8."""
    print_header("Python Version Check")
    version = sys.version_info
    print(f"Python: {version.major}.{version.minor}.{version.micro}")
    
    if version.major < 3 or (version.major == 3 and version.minor < 8):
        print("❌ Python 3.8+ required!")
        sys.exit(1)
    print("✅ Python version OK\n")


def install_dependencies():
    """Install requirements from requirements.txt."""
    print_header("Installing Dependencies")
    
    requirements_file = Path("requirements.txt")
    if not requirements_file.exists():
        print("❌ requirements.txt not found!")
        sys.exit(1)
    
    print("Installing packages from requirements.txt...")
    result = subprocess.run(
        [sys.executable, "-m", "pip", "install", "-r", "requirements.txt"],
        capture_output=False
    )
    
    if result.returncode != 0:
        print("❌ Installation failed!")
        sys.exit(1)
    
    print("✅ Dependencies installed\n")


def run_health_checks():
    """Verify key modules can be imported."""
    print_header("Health Checks")
    
    modules = [
        "streamlit",
        "plotly",
        "numpy",
        "pandas",
        "bioinformatics",
        "config",
        "export_utils",
    ]
    
    for module in modules:
        try:
            __import__(module)
            print(f"✅ {module}")
        except ImportError as e:
            print(f"❌ {module} — {e}")
            sys.exit(1)
    
    print("\n✅ All health checks passed!\n")


def check_database():
    """Verify gene database exists."""
    print_header("Database Check")
    
    db_path = Path("genes_database.json")
    if not db_path.exists():
        print(f"❌ {db_path} not found!")
        sys.exit(1)
    
    import json
    try:
        with open(db_path) as f:
            db = json.load(f)
        print(f"✅ Database loaded: {len(db)} genes")
        for gene_name in list(db.keys())[:3]:
            print(f"   - {gene_name}")
        print()
    except json.JSONDecodeError:
        print(f"❌ Error parsing {db_path}")
        sys.exit(1)


def launch_app():
    """Launch the Streamlit app."""
    print_header("Launching Application")
    
    print("Starting Streamlit app on http://localhost:8501")
    print("Press Ctrl+C to stop\n")
    
    subprocess.run(["streamlit", "run", "app.py"])


def main():
    """Main setup workflow."""
    print("\n")
    print("╔" + "=" * 68 + "╗")
    print("║" + " " * 68 + "║")
    print("║" + "  🧬 AI-Powered Plant Gene Analyzer - Setup & Launch".center(68) + "║")
    print("║" + "  v2.0 (May 2026)".center(68) + "║")
    print("║" + " " * 68 + "║")
    print("╚" + "=" * 68 + "╝\n")
    
    try:
        check_python_version()
        install_dependencies()
        run_health_checks()
        check_database()
        launch_app()
    
    except KeyboardInterrupt:
        print("\n\n❌ Setup interrupted by user")
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ Setup failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
