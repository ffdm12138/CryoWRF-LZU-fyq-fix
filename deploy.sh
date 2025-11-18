#!/bin/bash
# Deploy CRYOWRF-WRF from GitLab CI
set -euo pipefail

usage() {
    cat << EOF
Usage: ./deploy.sh [commit-sha] [-d DIR] [-l [N]] [-s]

Deploy CRYOWRF-WRF build artifacts from GitLab CI.

Arguments:
  commit-sha          Git commit hash (default: current HEAD)
  -d, --base-dir DIR  Installation directory (default: ./cryowrf-builds)
  -l, --list [N]      List recent pipelines with artifacts (default: 5)
  -s, --show          Show locally downloaded versions
  -h, --help          Show this help

Examples:
  ./deploy.sh                          # Deploy current commit
  ./deploy.sh abc1234                  # Deploy specific commit
  ./deploy.sh -d /scratch/wrf          # Custom location
  ./deploy.sh -l                       # Show last 5 pipelines with artifacts
  ./deploy.sh -l 10                    # Show last 10 pipelines with artifacts
  ./deploy.sh -s                       # Show downloaded versions

After deployment:
  cd cryowrf-builds/versions/<commit>/build
  export LD_LIBRARY_PATH=./lib:\$LD_LIBRARY_PATH
  export PATH=./bin:\$PATH
  wrf.exe
EOF
    exit 0
}

GITLAB_URL="https://gitlabext.wsl.ch"
GITLAB_WEB_URL="https://git.wsl.ch"
PROJECT_PATH="atmospheric-models/CRYOWRF-WRF"
PROJECT_ID="atmospheric-models%2FCRYOWRF-WRF"
BASE_DIR="${PWD}/cryowrf-builds"
COMMIT_SHA=""
LIST_COUNT=""

show_local() {
    echo "Downloaded versions:"
    echo ""

    if [ ! -d "${BASE_DIR}/versions" ]; then
        echo "No versions downloaded yet"
        echo ""
        echo "Use: ./deploy.sh <commit-sha> to download"
        exit 0
    fi

    # Header
    printf "%-8s  %-10s  %-10s  %s\n" \
        "COMMIT" "DATE" "SIZE" "PATH"
    printf "%-8s  %-10s  %-10s  %s\n" \
        "--------" "----------" "----------" "----"

    # List versions sorted by modification time (newest first)
    for version_dir in $(ls -t "${BASE_DIR}/versions" 2>/dev/null); do
        commit_short="$version_dir"
        version_path="${BASE_DIR}/versions/${version_dir}"

        # Get directory modification time and size
        mod_date=$(stat -c %y "$version_path" 2>/dev/null | cut -d' ' -f1)
        size=$(du -sh "$version_path" 2>/dev/null | cut -f1)

        printf "%-8s  %-10s  %-10s  %s\n" \
            "$commit_short" "$mod_date" "$size" "$version_path"
    done

    echo ""
    echo "To use a version:"
    echo "  cd ${BASE_DIR}/versions/<commit>/build"
    echo "  export LD_LIBRARY_PATH=./lib:\$LD_LIBRARY_PATH"
    echo "  export PATH=./bin:\$PATH"
    exit 0
}

list_pipelines() {
    local count="${1:-5}"
    echo "Available releases and recent builds:"
    echo ""

    # First, list published releases (from package registry)
    echo "=== Published Releases (permanent) ==="
    PACKAGES=$(curl -s "${GITLAB_URL}/api/v4/projects/${PROJECT_ID}/packages?package_type=generic&package_name=cryowrf-wrf&per_page=100&sort=desc")

    if [ -n "$PACKAGES" ] && [ "$PACKAGES" != "[]" ]; then
        printf "%-8s  %-10s  %-20s  %-11s  %s\n" \
            "TAG" "DATE" "NAME" "TYPE" "SIZE"
        printf "%-8s  %-10s  %-20s  %-11s  %s\n" \
            "--------" "----------" "--------------------" "-----------" "----"

        local shown_releases=0
        echo "$PACKAGES" | jq -r '.[] | "\(.version)|\(.created_at[0:10])|\(.name)"' | while IFS='|' read -r version date name; do
            [ $shown_releases -ge $count ] && break
            printf "%-8s  %-10s  %-20s  %-11s  %s\n" \
                "$version" "$date" "$name" "release" "-"
            shown_releases=$((shown_releases + 1))
        done
        echo ""
    else
        echo "No published releases yet"
        echo ""
    fi

    echo "=== Recent Builds (3-day artifacts) ==="
    # Get recent successful pipelines (more than requested to account for duplicates)
    PIPELINES=$(curl -s "${GITLAB_URL}/api/v4/projects/${PROJECT_ID}/pipelines?status=success&per_page=50")

    # Check if we got data
    if [ -z "$PIPELINES" ] || [ "$PIPELINES" = "[]" ]; then
        echo "No pipelines found"
        exit 1
    fi

    # Header
    printf "%-8s  %-10s  %-20s  %-11s  %s\n" \
        "COMMIT" "DATE" "BRANCH" "ARTIFACTS" "MESSAGE"
    printf "%-8s  %-10s  %-20s  %-11s  %s\n" \
        "--------" "----------" "--------------------" "-----------" "-------"

    declare -A seen_commits
    local shown=0

    while IFS='|' read -r pid sha ref created; do
        # Skip if we've already shown this commit (we only care about most recent pipeline)
        [ -n "${seen_commits[$sha]:-}" ] && continue
        seen_commits[$sha]=1

        [ $shown -ge $count ] && break

        # Check if this pipeline has build-wrf artifacts
        JOBS=$(curl -s "${GITLAB_URL}/api/v4/projects/${PROJECT_ID}/pipelines/${pid}/jobs")
        JOB_ID=$(echo "$JOBS" | jq -r '.[] | select(.name == "build-wrf" and .status == "success") | .id' | head -1)

        if [ -n "$JOB_ID" ]; then
            HAS_FILE=$(echo "$JOBS" | jq -r ".[] | select(.id == ${JOB_ID}) | .artifacts_file")

            # Check if artifacts actually exist (artifacts_file will be null if expired/unavailable)
            if [ "$HAS_FILE" = "null" ] || [ -z "$HAS_FILE" ]; then
                ARTIFACT_STATUS="✗ expired"
            else
                ARTIFACT_STATUS="✓"
            fi
        else
            ARTIFACT_STATUS="✗ no job"
        fi

        # Get commit info
        COMMIT_INFO=$(curl -s "${GITLAB_URL}/api/v4/projects/${PROJECT_ID}/repository/commits/${sha}")
        COMMIT_DATE=$(echo "$COMMIT_INFO" | jq -r '.created_at' | cut -d'T' -f1)
        MESSAGE=$(echo "$COMMIT_INFO" | jq -r '.title')

        # Truncate branch name if too long
        if [ ${#ref} -gt 20 ]; then
            ref="${ref:0:17}..."
        fi

        # Truncate message if too long
        if [ ${#MESSAGE} -gt 60 ]; then
            MESSAGE="${MESSAGE:0:57}..."
        fi

        printf "%-8s  %-10s  %-20s  %-11s  %s\n" \
            "${sha:0:7}" "$COMMIT_DATE" "$ref" "$ARTIFACT_STATUS" "$MESSAGE"

        shown=$((shown + 1))

    done < <(echo "$PIPELINES" | jq -r '.[] | "\(.id)|\(.sha)|\(.ref)|\(.created_at)"')

    echo ""
    echo "Deploy with:"
    echo "  ./deploy.sh <tag>         # For releases (e.g., v1.0.0)"
    echo "  ./deploy.sh <commit-sha>  # For recent builds"
    exit 0
}

# Parse args
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help) usage ;;
        -s|--show) show_local ;;
        -l|--list)
            shift
            if [[ $# -gt 0 ]] && [[ $1 =~ ^[0-9]+$ ]]; then
                list_pipelines "$1"
                shift
            else
                list_pipelines 5
            fi
            ;;
        -d|--base-dir) BASE_DIR="$2"; shift 2 ;;
        *) COMMIT_SHA="$1"; shift ;;
    esac
done

# Check prerequisites
command -v jq >/dev/null || { echo "Error: jq not found"; exit 1; }
command -v curl >/dev/null || { echo "Error: curl not found"; exit 1; }

# Get commit
if [ -z "$COMMIT_SHA" ]; then
    git rev-parse --git-dir >/dev/null 2>&1 || { echo "Error: not a git repo and no commit specified"; exit 1; }
    COMMIT_SHA=$(git rev-parse HEAD)
else
    # If a short SHA was provided, try to expand it via GitLab API
    if [ ${#COMMIT_SHA} -lt 40 ]; then
        echo "Looking up full commit SHA..."
        FULL_SHA=$(curl -s "${GITLAB_URL}/api/v4/projects/${PROJECT_ID}/repository/commits/${COMMIT_SHA}" | jq -r '.id // empty')
        if [ -n "$FULL_SHA" ]; then
            COMMIT_SHA="$FULL_SHA"
        fi
    fi
fi
COMMIT_SHORT=$(echo "$COMMIT_SHA" | cut -c1-7)

echo "Deploying CRYOWRF-WRF commit $COMMIT_SHORT to $BASE_DIR"

# Check if this is a tagged release in package registry
TAG=$(curl -s "${GITLAB_URL}/api/v4/projects/${PROJECT_ID}/repository/commits/${COMMIT_SHA}" | jq -r '.last_pipeline.ref // empty')
PACKAGE_URL=""

if [[ "$TAG" =~ ^v[0-9] ]]; then
    # Check if package exists for this tag
    PACKAGE=$(curl -s "${GITLAB_URL}/api/v4/projects/${PROJECT_ID}/packages?package_type=generic&package_name=cryowrf-wrf" | jq -r ".[] | select(.version == \"$TAG\")")
    if [ -n "$PACKAGE" ]; then
        echo "Found published release: $TAG"
        PACKAGE_URL="${GITLAB_URL}/api/v4/projects/${PROJECT_ID}/packages/generic/cryowrf-wrf/${TAG}/cryowrf-wrf-${TAG}.tar.gz"
    fi
fi

VERSION_DIR="${BASE_DIR}/versions/${COMMIT_SHORT}"
if [ -d "$VERSION_DIR" ]; then
    echo "Version $COMMIT_SHORT already exists"
else
    if [ -n "$PACKAGE_URL" ]; then
        # Download from package registry (permanent release)
        echo "Downloading release from package registry..."
        TMPFILE=$(mktemp)
        HTTP_CODE=$(curl -sL -w "%{http_code}" -o "$TMPFILE" "$PACKAGE_URL")

        if [ "$HTTP_CODE" != "200" ]; then
            rm -f "$TMPFILE"
            echo "Error: Failed to download package (HTTP $HTTP_CODE)"
            exit 1
        fi

        echo "Extracting release..."
        mkdir -p "$VERSION_DIR/build"
        if ! tar xzf "$TMPFILE" -C "$VERSION_DIR/build"; then
            rm -rf "$VERSION_DIR" "$TMPFILE"
            echo "Error: Failed to extract package"
            exit 1
        fi
        rm "$TMPFILE"
    else
        # Download from pipeline artifacts (3-day expiration)
        echo "Downloading from pipeline artifacts..."

        # Find pipeline
        PIPELINE=$(curl -s "${GITLAB_URL}/api/v4/projects/${PROJECT_ID}/pipelines?sha=${COMMIT_SHA}&per_page=1")
        PIPELINE_ID=$(echo "$PIPELINE" | jq -r '.[0].id // empty')
        STATUS=$(echo "$PIPELINE" | jq -r '.[0].status // empty')

        [ -z "$PIPELINE_ID" ] && { echo "Error: No pipeline found"; exit 1; }
        [ "$STATUS" != "success" ] && { echo "Error: Pipeline $STATUS (not success)"; exit 1; }

        echo "Found pipeline #$PIPELINE_ID"

        # Find build job
        JOBS=$(curl -s "${GITLAB_URL}/api/v4/projects/${PROJECT_ID}/pipelines/${PIPELINE_ID}/jobs")
        JOB_ID=$(echo "$JOBS" | jq -r '.[] | select(.name == "build-wrf" and .status == "success") | .id' | head -1)

        [ -z "$JOB_ID" ] && { echo "Error: No build-wrf job found"; exit 1; }

        EXPIRED=$(echo "$JOBS" | jq -r ".[] | select(.id == ${JOB_ID}) | .artifacts_expired")
        if [ "$EXPIRED" = "true" ]; then
            echo "Error: Artifacts expired (>3 days old)"
            echo ""
            echo "Re-run the 'build-wrf' job to rebuild artifacts:"
            echo "  ${GITLAB_WEB_URL}/${PROJECT_PATH}/-/pipelines/${PIPELINE_ID}"
            exit 1
        fi

        TMPFILE=$(mktemp)
        HTTP_CODE=$(curl -sL -w "%{http_code}" -o "$TMPFILE" \
            "${GITLAB_URL}/api/v4/projects/${PROJECT_ID}/jobs/${JOB_ID}/artifacts")

        if [ "$HTTP_CODE" != "200" ]; then
            rm -f "$TMPFILE"
            echo "Error: Failed to download artifacts (HTTP $HTTP_CODE)"
            echo ""
            echo "Re-run the 'build-wrf' job to rebuild artifacts:"
            echo "  ${GITLAB_WEB_URL}/${PROJECT_PATH}/-/pipelines/${PIPELINE_ID}"
            exit 1
        fi

        echo "Extracting..."
        mkdir -p "$VERSION_DIR"
        if ! unzip -q "$TMPFILE" -d "$VERSION_DIR"; then
            rm -rf "$VERSION_DIR" "$TMPFILE"
            echo "Error: Failed to extract artifacts"
            echo ""
            echo "Re-run the 'build-wrf' job to rebuild artifacts:"
            echo "  ${GITLAB_WEB_URL}/${PROJECT_PATH}/-/pipelines/${PIPELINE_ID}"
            exit 1
        fi
        rm "$TMPFILE"
    fi
fi

# Create env.sh with absolute paths (shell-agnostic)
cat > "$VERSION_DIR/env.sh" << EOF
#!/bin/bash
# CRYOWRF-WRF Environment Setup
# Source this file to set up paths for WRF executables and libraries

export CRYOWRF_VERSION="$COMMIT_SHORT"
export CRYOWRF_ROOT="$VERSION_DIR"
export LD_LIBRARY_PATH="$VERSION_DIR/build/lib:\${LD_LIBRARY_PATH}"
export PATH="$VERSION_DIR/build/bin:\${PATH}"

echo "CRYOWRF-WRF environment loaded (version: \${CRYOWRF_VERSION})"
echo "  Executables: wrf.exe, real.exe, ndown.exe, tc.exe"
echo "  Runtime files: \${CRYOWRF_ROOT}/build/run/"
EOF
chmod +x "$VERSION_DIR/env.sh"

echo ""
echo "Done! Deployed to $VERSION_DIR"
echo ""
echo "To use:"
echo "  source $VERSION_DIR/env.sh"
echo "  cd <your-working-directory>  # Directory with namelist.input"
echo "  real.exe && wrf.exe"
