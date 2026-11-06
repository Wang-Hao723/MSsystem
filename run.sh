#!/bin/bash
# MSsystem V16 — one logical UI/API endpoint with safe process ownership.
#
# Usage:
#   bash run.sh                         start V16 on 127.0.0.1:8516
#   bash run.sh --port 8600             explicit alternate port
#   bash run.sh --foreground            keep server attached for debugging
#   bash run.sh --status                inspect the managed V16 process
#   bash run.sh --stop                  stop only this checkout's process
#   bash run.sh --small --no-open       use the small DB and do not open browser
#   bash run.sh --streamlit             explicit legacy UI on port 8526

set -u

ACTION="start"
MODE="fastapi"
HOST="${MSSYSTEM_HOST:-127.0.0.1}"
PORT="${MSSYSTEM_PORT:-8516}"
NO_OPEN=0
FOREGROUND=0
START_WAIT_SECONDS="${MSSYSTEM_START_WAIT_SECONDS:-1800}"
PROJECT_ROOT="$(cd "$(dirname "$0")" && pwd)"

resolve_python() {
    if [[ -n "${MSSYSTEM_PYTHON:-}" ]]; then
        printf '%s\n' "$MSSYSTEM_PYTHON"
        return 0
    fi
    if [[ -n "${CONDA_PREFIX:-}" ]] && [[ -x "${CONDA_PREFIX}/bin/python" ]]; then
        printf '%s\n' "${CONDA_PREFIX}/bin/python"
        return 0
    fi
    # This workstation's validated MSsystem runtime lives outside the shell's
    # default PATH.  Prefer it over an unrelated framework python3 that may not
    # have FastAPI, uvicorn, RDKit, or the mass-spectrometry dependencies.
    if [[ -x "/opt/anaconda3/bin/python" ]]; then
        printf '%s\n' "/opt/anaconda3/bin/python"
        return 0
    fi
    if command -v python3 >/dev/null 2>&1; then
        command -v python3
        return 0
    fi
    if command -v python >/dev/null 2>&1; then
        command -v python
        return 0
    fi
    return 1
}

if ! PYTHON_BIN="$(resolve_python)"; then
    echo "No Python interpreter found. Activate the MSsystem environment or set MSSYSTEM_PYTHON." >&2
    exit 2
fi
RUNTIME_STATE_DIR="${PROJECT_ROOT}/output/runtime"
PID_FILE="${RUNTIME_STATE_DIR}/mssystem_v16.pid"
ENDPOINT_FILE="${RUNTIME_STATE_DIR}/mssystem_v16.endpoint"

usage() {
    sed -n '2,12p' "$0"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --status)
            ACTION="status"
            shift
            ;;
        --stop)
            ACTION="stop"
            shift
            ;;
        --streamlit)
            MODE="streamlit"
            PORT="8526"
            shift
            ;;
        --small)
            export MSSYSTEM_DB_SMALL=1
            shift
            ;;
        --foreground)
            FOREGROUND=1
            shift
            ;;
        --no-open)
            NO_OPEN=1
            shift
            ;;
        --port)
            [[ $# -ge 2 ]] || { echo "Error: --port requires a value" >&2; exit 2; }
            PORT="$2"
            shift 2
            ;;
        --host)
            [[ $# -ge 2 ]] || { echo "Error: --host requires a value" >&2; exit 2; }
            HOST="$2"
            shift 2
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            echo "Error: unknown option $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

case "$PORT" in
    ''|*[!0-9]*)
        echo "Error: port must be an integer" >&2
        exit 2
        ;;
esac
if (( PORT < 1 || PORT > 65535 )); then
    echo "Error: port must be between 1 and 65535" >&2
    exit 2
fi

# Source-bound full-library verification can exceed three minutes on the
# external disk. Keep a finite, configurable deadline for the actual process.
if [[ "$ACTION" = "start" ]]; then
    if [[ ! "$START_WAIT_SECONDS" =~ ^[0-9]{1,4}$ ]] || (( 10#$START_WAIT_SECONDS < 1 || 10#$START_WAIT_SECONDS > 3600 )); then
        echo "Error: MSSYSTEM_START_WAIT_SECONDS must be an integer from 1 to 3600" >&2
        exit 2
    fi
    START_WAIT_SECONDS=$((10#$START_WAIT_SECONDS))
fi

managed_pid() {
    [[ -f "$PID_FILE" ]] || return 1
    local pid
    pid="$(tr -cd '0-9' < "$PID_FILE")"
    [[ -n "$pid" ]] || return 1
    kill -0 "$pid" 2>/dev/null || return 1
    printf '%s\n' "$pid"
}

pid_cwd() {
    local pid="$1"
    if [[ -e "/proc/${pid}/cwd" ]]; then
        readlink "/proc/${pid}/cwd" 2>/dev/null
        return
    fi
    if command -v lsof >/dev/null 2>&1; then
        lsof -a -p "$pid" -d cwd -Fn 2>/dev/null | sed -n 's/^n//p' | head -n 1
    fi
}

pid_is_ours() {
    local pid="$1"
    local cwd command
    cwd="$(pid_cwd "$pid")"
    command="$(ps -p "$pid" -o command= 2>/dev/null)"
    [[ "$cwd" = "$PROJECT_ROOT" ]] || return 1
    [[ "$command" == *"web.server:app"* || "$command" == *"web/app.py"* ]]
}

show_status() {
    local pid status_host status_port status_mode
    status_host="$HOST"
    status_port="$PORT"
    if [[ -f "$ENDPOINT_FILE" ]]; then
        read -r status_host status_port < "$ENDPOINT_FILE"
    fi
    if pid="$(managed_pid)" && pid_is_ours "$pid"; then
        echo "PID: $pid"
        echo "Project: $PROJECT_ROOT"
        echo "Endpoint: http://${status_host}:${status_port}"
        status_mode="fastapi"
        if [[ "$(ps -p "$pid" -o command= 2>/dev/null)" == *"streamlit"* ]]; then
            status_mode="streamlit"
        fi
        if "$PYTHON_BIN" "${PROJECT_ROOT}/scripts/service_launcher.py" --check \
            --project-root "$PROJECT_ROOT" --host "$status_host" --port "$status_port" --mode "$status_mode"; then
            echo "MSsystem V16 is running and healthy"
            return 0
        fi
        echo "Process exists but the endpoint is not ready; inspect output/runtime/mssystem_v16_service.log" >&2
        return 4
    fi
    echo "MSsystem V16 is not running from this checkout"
    return 1
}

mkdir -p "$RUNTIME_STATE_DIR"

if [[ "$ACTION" = "status" ]]; then
    show_status
    exit $?
fi

if [[ "$ACTION" = "stop" ]]; then
    if ! pid="$(managed_pid)"; then
        echo "MSsystem V16 has no live managed PID"
        exit 0
    fi
    if ! pid_is_ours "$pid"; then
        echo "Refusing to stop PID $pid: it is not owned by $PROJECT_ROOT" >&2
        exit 3
    fi
    kill "$pid"
    echo "Graceful stop requested for MSsystem V16 PID $pid"
    # Uvicorn may close its listener before Python releases a multi-GB library
    # and executor threads.  Do not report a completed stop while the managed
    # PID is still alive: a subsequent start would otherwise be refused as
    # "already running" even though the UI is unreachable.
    for _ in $(seq 1 150); do
        if ! kill -0 "$pid" 2>/dev/null; then
            echo "MSsystem V16 stopped"
            exit 0
        fi
        sleep 0.1
    done
    if ! pid_is_ours "$pid"; then
        echo "Refusing forced stop: PID $pid no longer belongs to $PROJECT_ROOT" >&2
        exit 3
    fi
    echo "Graceful stop timed out; forcing only verified V16 PID $pid to exit" >&2
    kill -KILL "$pid"
    for _ in $(seq 1 50); do
        if ! kill -0 "$pid" 2>/dev/null; then
            echo "MSsystem V16 stopped after forced process cleanup"
            exit 0
        fi
        sleep 0.1
    done
    echo "MSsystem V16 PID $pid did not exit" >&2
    exit 4
fi

if pid="$(managed_pid)"; then
    if pid_is_ours "$pid"; then
        echo "MSsystem V16 is already running as PID $pid"
        show_status
        exit $?
    fi
    mv "$PID_FILE" "${PID_FILE}.stale.$(date +%Y%m%d%H%M%S)"
fi

if command -v lsof >/dev/null 2>&1; then
    occupied="$(lsof -tiTCP:"${PORT}" -sTCP:LISTEN 2>/dev/null | head -n 1)"
    if [[ -n "$occupied" ]]; then
        echo "Port ${PORT} is already occupied by PID ${occupied}." >&2
        ps -p "$occupied" -o pid=,command= 2>/dev/null >&2 || true
        echo "Choose --port <other>, or stop that program yourself. V16 will not kill it." >&2
        exit 2
    fi
fi

if [[ ! -x "$PYTHON_BIN" ]]; then
    echo "Configured Python is not executable: $PYTHON_BIN" >&2
    exit 2
fi
if ! "$PYTHON_BIN" -c "import fastapi, uvicorn" >/dev/null 2>&1; then
    echo "FastAPI/uvicorn are unavailable in $PYTHON_BIN; dependencies were not auto-installed." >&2
    exit 2
fi

export MSSYSTEM_HOST="$HOST"
export MSSYSTEM_PORT="$PORT"
export MSSYSTEM_RUNTIME_ROOT="${MSSYSTEM_RUNTIME_ROOT:-${PROJECT_ROOT}/output}"
# Converted vendor RAW files are reusable across V16 restarts and can be large.
# Prefer the user-mounted high-capacity disk when it is present; an explicit
# environment setting always wins, and the runtime directory is the safe
# fallback on another workstation.
if [[ -z "${MSSYSTEM_VENDOR_CONVERSION_CACHE_ROOT:-}" ]]; then
    if [[ -d "/Volumes/LaCie" ]] && [[ -w "/Volumes/LaCie" ]]; then
        export MSSYSTEM_VENDOR_CONVERSION_CACHE_ROOT="/Volumes/LaCie/Mssystem_V16/vendor_conversion_cache"
    else
        export MSSYSTEM_VENDOR_CONVERSION_CACHE_ROOT="${MSSYSTEM_RUNTIME_ROOT}/vendor_conversion_cache"
    fi
fi
export MPLCONFIGDIR="${MPLCONFIGDIR:-${PROJECT_ROOT}/output/.matplotlib}"
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-${PROJECT_ROOT}/output/.numba}"
export PYTHONPYCACHEPREFIX="${PYTHONPYCACHEPREFIX:-/tmp/mssystem_v16_pycache}"
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR" "$MSSYSTEM_VENDOR_CONVERSION_CACHE_ROOT"

"$PYTHON_BIN" "${PROJECT_ROOT}/scripts/storage_preflight.py" || exit $?
"$PYTHON_BIN" "${PROJECT_ROOT}/scripts/generate_runtime_manifest.py" \
    --source-only --manifest "${PROJECT_ROOT}/config/runtime_manifest.json" || exit $?

echo "═══════════════════════════════════════════════════════════"
echo "  MSsystem V16 / V15-Flagship + Dual-Polarity Shadow Evidence"
echo "  Project: $PROJECT_ROOT"
echo "  UI + API: http://${HOST}:${PORT}"
echo "  API docs: http://${HOST}:${PORT}/docs"
echo "  Process control: bash run.sh --status | --stop"
echo "═══════════════════════════════════════════════════════════"

if [[ "$FOREGROUND" -eq 0 ]]; then
    "$PYTHON_BIN" "${PROJECT_ROOT}/scripts/service_launcher.py" \
        --project-root "$PROJECT_ROOT" --python "$PYTHON_BIN" \
        --host "$HOST" --port "$PORT" --mode "$MODE" --wait-seconds "$START_WAIT_SECONDS"
    start_status=$?
    if [[ "$start_status" -eq 0 ]] && [[ "$NO_OPEN" -eq 0 ]]; then
        if command -v open >/dev/null 2>&1; then
            open "http://${HOST}:${PORT}" >/dev/null 2>&1 || true
        elif command -v xdg-open >/dev/null 2>&1; then
            xdg-open "http://${HOST}:${PORT}" >/dev/null 2>&1 || true
        fi
    fi
    exit "$start_status"
fi

if [[ "$NO_OPEN" -eq 0 ]] && command -v open >/dev/null 2>&1; then
    (sleep 2; open "http://${HOST}:${PORT}") >/dev/null 2>&1 &
elif [[ "$NO_OPEN" -eq 0 ]] && command -v xdg-open >/dev/null 2>&1; then
    (sleep 2; xdg-open "http://${HOST}:${PORT}") >/dev/null 2>&1 &
fi

if [[ "$MODE" = "streamlit" ]]; then
    "$PYTHON_BIN" -m streamlit run web/app.py \
        --server.port "$PORT" --server.address "$HOST" &
else
    "$PYTHON_BIN" -m uvicorn web.server:app \
        --host "$HOST" --port "$PORT" &
fi

SERVER_PID=$!
printf '%s\n' "$SERVER_PID" > "$PID_FILE"
printf '%s %s\n' "$HOST" "$PORT" > "$ENDPOINT_FILE"

cleanup() {
    if [[ -f "$PID_FILE" ]] && [[ "$(tr -cd '0-9' < "$PID_FILE")" = "$SERVER_PID" ]]; then
        suffix="$(date +%Y%m%d%H%M%S)"
        mv "$PID_FILE" "${PID_FILE}.stopped.${suffix}"
        if [[ -f "$ENDPOINT_FILE" ]]; then
            mv "$ENDPOINT_FILE" "${ENDPOINT_FILE}.stopped.${suffix}"
        fi
    fi
}
trap cleanup EXIT
wait "$SERVER_PID"
