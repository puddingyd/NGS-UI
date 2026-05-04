"""Bootstrap CLI: `python -m app create-user [USERNAME]`, `list-users`.

Run from the backend/ parent dir with PYTHONPATH=backend.
"""
from __future__ import annotations

import getpass
import sys

from .services import users


def _create_user(argv: list[str]) -> int:
    if argv:
        username = argv[0]
    else:
        username = input("username: ").strip()
    password = getpass.getpass("password: ")
    confirm  = getpass.getpass("confirm:  ")
    if password != confirm:
        print("ERROR: passwords don't match", file=sys.stderr)
        return 2
    try:
        u = users.create_user(username, password)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr); return 2
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr); return 1
    print(f"created user id={u['id']}  username={u['username']}")
    return 0


def _list_users() -> int:
    rows = users.list_users()
    if not rows:
        print("(no users yet — run create-user)")
        return 0
    for r in rows:
        print(f"  id={r['id']:<3}  {r['username']:<20}  created {r['created_at']}")
    return 0


def main(argv: list[str]) -> int:
    if not argv or argv[0] in ("-h", "--help"):
        print(
            "usage: python -m app COMMAND [...]\n"
            "  create-user [USERNAME]   create a new account\n"
            "  list-users               list existing accounts\n",
            file=sys.stderr,
        )
        return 0
    cmd, rest = argv[0], argv[1:]
    if cmd == "create-user":
        return _create_user(rest)
    if cmd == "list-users":
        return _list_users()
    print(f"unknown command: {cmd}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
