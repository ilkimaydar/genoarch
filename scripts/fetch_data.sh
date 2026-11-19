#!/usr/bin/env sh
# Download the pinned HPO release files (and optionally Phenopacket Store) and verify sha256.
# Usage: scripts/fetch_data.sh [--with-store]
set -eu
HPO_TAG="${HPO_TAG:-v2026-09-01}"
STORE_TAG="${STORE_TAG:-0.1.27}"
DEST="${DEST:-data}"
mkdir -p "$DEST/hpo"
base="https://github.com/obophenotype/human-phenotype-ontology/releases/download/$HPO_TAG"
for f in hp.json hp.obo phenotype.hpoa genes_to_disease.txt genes_to_phenotype.txt; do
  [ -s "$DEST/hpo/$f" ] || curl -fsSL -o "$DEST/hpo/$f" "$base/$f"
done
if [ "$HPO_TAG" = "v2026-09-01" ]; then
  (cd "$DEST" && sha256sum -c - <<'SUMS'
a7b3a012e7b4007a35a7cf8da35f2373b54cc16907f80d13d62470d31f833501  hpo/hp.json
93dace952fcb3ec4728818857f6ba76bc2d5312f4d83266519b8694b8e798f22  hpo/hp.obo
e89aa39c8f97bf5a52c5d160f9250a603681d630ade8ec2679d84ac5aece1f72  hpo/phenotype.hpoa
507a17bff9c49e6329fbd88b1f91734fa7e9fdea5c06304044c0892eb6ab248c  hpo/genes_to_phenotype.txt
7937a1b92779dc8614e39ace623b443c7c71dd15e4de11d5e52b14e331ab5978  hpo/genes_to_disease.txt
SUMS
  )
fi
if [ "${1:-}" = "--with-store" ]; then
  mkdir -p "$DEST/phenopacket-store/$STORE_TAG"
  zip="$DEST/phenopacket-store/$STORE_TAG/all_phenopackets.zip"
  [ -s "$zip" ] || curl -fsSL -o "$zip" "https://github.com/monarch-initiative/phenopacket-store/releases/download/$STORE_TAG/all_phenopackets.zip"
  if [ "$STORE_TAG" = "0.1.27" ]; then
    echo "d0c70005bb09b87035087516252b6b539fbc0415be7c664b73c9cde7eedf205a  $zip" | sha256sum -c -
  fi
  python3 -c "import zipfile,sys; zipfile.ZipFile(sys.argv[1]).extractall(sys.argv[2])" "$zip" "$DEST/phenopacket-store/$STORE_TAG/extracted"
fi
echo "data ready in $DEST"
