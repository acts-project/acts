#!/bin/bash
set -e
set -u

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
export GIT_CLIFF_CONFIG=$SCRIPT_DIR/cliff.toml

# helper function to selectively print and run commands without a subshell
function run() {
    set -x
    "$@"
    # save exit code
    { rec=$?; } 2> /dev/null
    { set +x;   } 2> /dev/null
    # restore exit code
    (exit $rec)
}

export run

if [ -z "${version:-}" ]; then
  version=$(git cliff --bumped-version)
fi
echo "Bumped version will be: $version"

zenodo=$(cat .zenodo.json)
echo "$zenodo" \
  | jq --arg version "$version" '.version = $version' \
  | jq --arg version "$version" '.title = "acts-project/acts: \($version)"' \
  > .zenodo.json
echo "- Updated .zenodo.json"

citation=$(cat CITATION.cff)
echo "$citation" \
  | sed "s/^version: .*/version: $version/" \
  > CITATION.cff
echo "- Updated CITATION.cff"

echo "$version" | sed 's/^v//g' > version_number
echo "- Updated version_bumber"

run git add .zenodo.json CITATION.cff version_number
run git commit -n -m"Bump version to $version"

RELEASE_TARGET=${RELEASE_TARGET:-$(git rev-parse HEAD)}

# Build and verify the source archive before anything is pushed, so that a
# problem with it leaves no bumped commit and no release behind. Building it
# locally rather than downloading it from GitHub is what makes that ordering
# possible, since there is nothing to download until the commit is pushed.
#
# This archive is not the one GitHub generates for the tag: it additionally
# carries the generated code in `prebuilt-codegen/`, which lets a build from it
# skip the code generators and the downloads they need.
archive="acts-${version}.tar.gz"
archive_prefix="acts-${version#v}"
run uv run --no-project "$SCRIPT_DIR/pregenerate_codegen.py" \
  --output prebuilt-codegen
run "$SCRIPT_DIR/make_source_archive.sh" \
  --output "$archive" \
  --prebuilt prebuilt-codegen \
  --revision "$RELEASE_TARGET" \
  --prefix "$archive_prefix"

# Check the archive the way a consumer sees it: unpacked, with nothing pointing
# at the generated code, and with any fall back to running a generator fatal.
rm -rf release-verify
mkdir release-verify
run tar xzf "$archive" -C release-verify
run cmake -S "release-verify/${archive_prefix}/CI/codegen_prebuilt_check" \
  -B release-verify/build
run cmake --build release-verify/build
rm -rf release-verify prebuilt-codegen

CI=${CI:-}
if [ -n "$CI" ]; then
  run git push
fi

run git cliff --tag "$version" --latest --unreleased -o release.md

set +e
! gh release view "$version" > /dev/null 2>&1
release_exists=$?
set -e
if [[ $release_exists == 1 ]]; then
  echo "Release $version exists"
  run gh release edit $version \
    --notes-file release.md \
    --target $RELEASE_TARGET \
    --draft
else
  echo "Release $version does not exist"
  run gh release create $version \
    --title "$version" \
    --notes-file release.md \
    --target $RELEASE_TARGET \
    --draft
fi

run gh release upload "${version}" "acts-${version}.tar.gz" \
  --clobber
