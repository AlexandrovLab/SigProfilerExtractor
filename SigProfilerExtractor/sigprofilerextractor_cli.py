#!/usr/bin/env python3

import sys
from SigProfilerExtractor.controllers import cli_controller


def main_function():
    commands = {
        "sigprofilerextractor": "Extract mutational signatures from input samples."
    }

    if len(sys.argv) < 2 or sys.argv[1].lower() not in commands:
        print_usage(commands)
        sys.exit(1)

    command = sys.argv[1].lower()
    args = sys.argv[2:]

    controller = cli_controller.CliController()

    if command == "sigprofilerextractor":
        controller.dispatch_sigProfilerExtractor(args)


def print_usage(commands):
    """Prints the usage message."""
    print("Usage: SigProfilerExtractor <command> [<args>]\n")
    print("Commands:")
    for cmd, desc in commands.items():
        print(f"  {cmd}: {desc}")


if __name__ == "__main__":
    main_function()
