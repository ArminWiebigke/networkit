#!/bin/bash
for FILE in $1*/*.pdf; do
  echo ${FILE}
  pdfcrop  --margins '0 0 0 0' "${FILE}" "${FILE}"
done
