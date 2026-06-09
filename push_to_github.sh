#!/bin/zsh
# ─────────────────────────────────────────────────────────────
# push_to_github.sh
# Creates the TARVIA-lab/HGSOC-virtual-screening repo and pushes.
# Requires a GitHub Personal Access Token with `repo` scope.
# Get one at: https://github.com/settings/tokens/new
# ─────────────────────────────────────────────────────────────

set -e

echo ""
echo "══════════════════════════════════════════════════════"
echo "  HGSOC Virtual Screening → GitHub Push Script"
echo "══════════════════════════════════════════════════════"
echo ""

# ── 1. Collect credentials ────────────────────────────────────
read "GH_USER?GitHub username (e.g. omarlujanoolazaba): "
read -s "GH_TOKEN?Personal Access Token (repo scope): "
echo ""

REPO_NAME="HGSOC-virtual-screening"
ORG="TARVIA-lab"
DESCRIPTION="AI-driven drug repurposing for RelB-dependent targets in HGSOC"

echo ""
echo "── Step 1/4  Creating repo ${ORG}/${REPO_NAME} via GitHub API..."
HTTP_CODE=$(curl -s -o /tmp/gh_api_response.json -w "%{http_code}" \
  -X POST \
  -H "Authorization: token ${GH_TOKEN}" \
  -H "Accept: application/vnd.github.v3+json" \
  "https://api.github.com/orgs/${ORG}/repos" \
  -d "{
    \"name\": \"${REPO_NAME}\",
    \"description\": \"${DESCRIPTION}\",
    \"private\": false,
    \"auto_init\": false
  }")

if [ "$HTTP_CODE" = "201" ]; then
  echo "   ✓ Repository created: https://github.com/${ORG}/${REPO_NAME}"
elif [ "$HTTP_CODE" = "422" ]; then
  echo "   ℹ  Repository already exists — continuing with push."
else
  echo "   ✗ GitHub API returned HTTP ${HTTP_CODE}:"
  cat /tmp/gh_api_response.json
  echo ""
  echo "   → Make sure your token has the 'repo' scope and belongs to TARVIA-lab."
  exit 1
fi

# ── 2. Set remote ─────────────────────────────────────────────
echo ""
echo "── Step 2/4  Configuring git remote..."
cd "$(dirname "$0")"

if git remote get-url origin &>/dev/null; then
  git remote set-url origin "https://${GH_USER}:${GH_TOKEN}@github.com/${ORG}/${REPO_NAME}.git"
  echo "   ✓ Remote updated."
else
  git remote add origin "https://${GH_USER}:${GH_TOKEN}@github.com/${ORG}/${REPO_NAME}.git"
  echo "   ✓ Remote added."
fi

# ── 3. Push ───────────────────────────────────────────────────
echo ""
echo "── Step 3/4  Pushing to GitHub..."
git push -u origin main
echo "   ✓ Push complete."

# ── 4. Open in browser ────────────────────────────────────────
echo ""
echo "── Step 4/4  Opening repository in browser..."
REPO_URL="https://github.com/${ORG}/${REPO_NAME}"
open "$REPO_URL"

echo ""
echo "══════════════════════════════════════════════════════"
echo "  Done!  ${REPO_URL}"
echo "══════════════════════════════════════════════════════"
echo ""
