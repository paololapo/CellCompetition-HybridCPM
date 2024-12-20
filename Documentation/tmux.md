# Using `tmux` for Persistent Terminal Sessions
`tmux` is a terminal multiplexer that allows you to run multiple terminal sessions within a single window. It is especially useful when working on remote virtual machines, as it allows you to keep processes running even after closing the terminal. Here's how to use `tmux` effectively:

### 1. Creating a New `tmux` Session
To start a new `tmux` session, use the following command:
```bash
tmux new -s <session_name>
```
Replace `<session_name>` with a name of your choice. This command creates a new session and attaches you to it.

### 2. Listing Existing `tmux` Sessions
If you have multiple `tmux` sessions running and want to see them, list the existing sessions with:
```bash
tmux ls
```
This command will display all active sessions, their names, and their IDs.

### 3. Attaching to an Existing `tmux` Session
To reattach to a `tmux` session that you've previously detached from, use:
```bash
tmux attach -t <session_name>
```
Replace `<session_name>` with the name of the session you want to reattach to. If you have only one session, you can simply use:
```bash
tmux attach
```

### 4. Detaching from a `tmux` Session
If you want to leave a `tmux` session running while closing the terminal or switching to another session, detach from it by pressing:
```
Ctrl + b, then d
```
This will leave the session running in the background.

### 5. Scrolling Through Terminal Output in a `tmux` Session
To scroll through the terminal output within a `tmux` session:
1. Enter copy mode by pressing:
   ```
   Ctrl + b, then [
   ```
2. Use the arrow keys or `PgUp`/`PgDn` to scroll through the output.
3. To exit copy mode, press:
   ```
   q
   ```
