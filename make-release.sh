#!/usr/bin/env bash
#
# make-release.sh
#
# Creates a release tarball pluto-<version>.tar.gz in the current directory,
# where <version> is the one set in AC_INIT in configure.ac. The tarball
# contains the sources at HEAD of this repository along with all submodules
# (at their recorded commits) and pre-generated autotools files so that a user
# only needs to run ./configure && make.
#
# Usage: ./make-release.sh [-f]
#   -f  Allow creating a release even if the working tree has uncommitted
#       changes to tracked files (they are NOT included in the release).

set -euo pipefail

force=0
if [ "${1:-}" = "-f" ]; then
  force=1
fi

srcdir=$(cd "$(dirname "$0")" && pwd)
cd "$srcdir"

version=$(sed -n 's/^AC_INIT(\[pluto\],[[:space:]]*\[\([^]]*\)\].*/\1/p' configure.ac)
if [ -z "$version" ]; then
  echo "Error: could not determine version from configure.ac" >&2
  exit 1
fi

name=pluto-$version
tarball=$srcdir/$name.tar.gz

if [ $force -eq 0 ] && ! git diff-index --quiet HEAD --ignore-submodules=untracked; then
  echo "Error: working tree has uncommitted changes; commit them or use -f." >&2
  exit 1
fi

# Submodules must be checked out at the commits recorded in HEAD.
if git submodule status --recursive | grep -q '^[-+U]'; then
  echo "Error: submodules are not initialized or not at their recorded commits:" >&2
  git submodule status --recursive | grep '^[-+U]' >&2
  echo "Run 'git submodule update --init'." >&2
  exit 1
fi

for tool in autoreconf automake libtoolize; do
  if ! command -v $tool > /dev/null; then
    echo "Error: $tool is required to create a release." >&2
    exit 1
  fi
done

workdir=$(mktemp -d)
trap 'rm -rf "$workdir"' EXIT
stagedir=$workdir/$name

echo "*** Exporting pluto $version sources ***"
mkdir -p "$stagedir"
git archive --format=tar HEAD | tar -x -C "$stagedir"

git submodule foreach --quiet --recursive 'echo "$displaypath"' | while read -r sm; do
  echo "*** Exporting submodule $sm ***"
  mkdir -p "$stagedir/$sm"
  (cd "$sm" && git archive --format=tar HEAD) | tar -x -C "$stagedir/$sm"
  # isl and pet derive their version from git; record it since the release
  # has no .git.
  (cd "$sm" && git describe --always) > "$stagedir/$sm/GIT_HEAD_ID"
done

echo "*** Generating autotools files ***"
(cd "$stagedir" && ./autogen.sh)

# Remove files not meant for distribution.
find "$stagedir" -name autom4te.cache -type d -prune -exec rm -rf {} +
find "$stagedir" \( -name .gitignore -o -name .gitmodules -o -name .gitattributes \) \
  -type f -delete
rm -rf "$stagedir/.github"

echo "*** Creating $name.tar.gz ***"
tar -C "$workdir" -czf "$tarball" --owner=0 --group=0 "$name"

echo "Created $tarball"
