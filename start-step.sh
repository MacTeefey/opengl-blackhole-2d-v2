#!/bin/bash

# start-step.sh - Copy starter files for a specific step
# Usage: ./start-step.sh <step-folder-name>
# Example: ./start-step.sh 01a-cmake-setup

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STEPS_DIR="$SCRIPT_DIR/steps"

# Check if step name was provided
if [ -z "$1" ]; then
    echo "Usage: ./start-step.sh <step-folder-name>"
    echo ""
    echo "Available steps:"
    ls -1 "$STEPS_DIR" | sort
    exit 1
fi

STEP_NAME="$1"
STEP_PATH="$STEPS_DIR/$STEP_NAME"

# Check if step exists
if [ ! -d "$STEP_PATH" ]; then
    echo "Error: Step '$STEP_NAME' not found."
    echo ""
    echo "Available steps:"
    ls -1 "$STEPS_DIR" | sort
    exit 1
fi

echo "Starting step: $STEP_NAME"
echo "================================"

# Check for uncommitted changes in files that would be overwritten
if [ -d ".git" ]; then
    MODIFIED_FILES=""
    shopt -s dotglob
    for file in "$STEP_PATH"/*; do
        if [ -f "$file" ]; then
            filename=$(basename "$file")
            if [ -f "$SCRIPT_DIR/$filename" ]; then
                # Check if file has uncommitted changes
                if git diff --name-only "$SCRIPT_DIR/$filename" 2>/dev/null | grep -q . || \
                   git diff --staged --name-only "$SCRIPT_DIR/$filename" 2>/dev/null | grep -q .; then
                    MODIFIED_FILES="$MODIFIED_FILES  - $filename\n"
                fi
            fi
        fi
    done
    shopt -u dotglob

    if [ -n "$MODIFIED_FILES" ]; then
        echo ""
        echo "Warning: The following files have uncommitted changes that will be overwritten:"
        echo -e "$MODIFIED_FILES"
        read -p "Continue and overwrite these files? [y/N] " -n 1 -r
        echo ""
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            echo "Aborted. No files were changed."
            exit 0
        fi
    fi
fi

# Copy all files from the step directory to current directory (including hidden files)
shopt -s dotglob
for file in "$STEP_PATH"/*; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        cp "$file" "$SCRIPT_DIR/$filename"
        echo "  Copied: $filename"
    fi
done
shopt -u dotglob

# Colors for output
GREEN='\033[0;32m'
NC='\033[0m' # No Color

# Enable git hooks if .githooks directory exists
if [ -d ".git" ] && [ -d ".githooks" ]; then
    CURRENT_HOOKS_PATH=$(git config --local core.hooksPath 2>/dev/null || echo "")
    if [ "$CURRENT_HOOKS_PATH" != ".githooks" ]; then
        git config --local core.hooksPath .githooks
        echo ""
        echo -e "${GREEN}✓${NC} Git hooks enabled (pre-push compilation check)"
    fi
fi

echo ""
echo "Step '$STEP_NAME' starter files are now in place."
echo ""
echo "To build and run locally:"
echo "  ./run.sh"
echo ""
echo "To build on HelloC++:"
echo "  git add ."
echo "  git commit -m \"Start $STEP_NAME\""
echo "  git push"
echo ""
echo "Stuck? Send a screenshot and description of your problem to"
echo "hellocppdotdev@gmail.com for help."
