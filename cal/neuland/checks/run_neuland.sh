#!/bin/bash
#
# run_neuland.sh -- lanza neulandCheck.C en paralelo sobre todos los ficheros
#                   de files.txt y fusiona las salidas con hadd.
#
# Uso:
#   ./run_neuland.sh [n_cores] [files.txt] [outdir] [merged.root]
#
# Ejemplos:
#   ./run_neuland.sh                 # 8 cores, files.txt, out_reduced/, neuland_reduced_all.root
#   ./run_neuland.sh 16              # 16 cores
#   ./run_neuland.sh 16 mis_files.txt /scratch/salidas todo.root
#
# Variables de entorno opcionales:
#   MACRO=neulandCheck.C   nombre/ruta de la macro
#   FORCE=1                reprocesa aunque el .root de salida ya exista
#   MERGE_ONLY=1           salta el procesado y solo hace el hadd
#   NO_MERGE=1             procesa pero no hace el hadd
#

set -uo pipefail

# ----------------------------- configuracion -----------------------------
NCORES=${1:-8}
FILELIST=${2:-files.txt}
OUTDIR=${3:-out_reduced}
MERGED=${4:-neuland_reduced_all.root}

MACRO=${MACRO:-neulandCheck.C}
LOGDIR="${OUTDIR}/logs"

# Nombre de la funcion = nombre de la macro sin extension ni ruta
MACRO_ABS=$(readlink -f "$MACRO" 2>/dev/null || echo "$MACRO")

# ------------------------------ comprobaciones ---------------------------
command -v root >/dev/null 2>&1 || { echo "ERROR: 'root' no esta en el PATH. Carga tu entorno de ROOT/R3BRoot."; exit 1; }
command -v hadd >/dev/null 2>&1 || { echo "ERROR: 'hadd' no esta en el PATH."; exit 1; }
[[ -f "$MACRO_ABS" ]] || { echo "ERROR: no encuentro la macro '$MACRO'."; exit 1; }
[[ -f "$FILELIST"  ]] || { echo "ERROR: no encuentro la lista '$FILELIST'."; exit 1; }

mkdir -p "$OUTDIR" "$LOGDIR"

# ------------------------- lectura de la lista ---------------------------
# Acepta lineas con o sin comillas, ignora vacias y comentarios (#).
mapfile -t INFILES < <(sed -e 's/\r$//' -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' \
                           -e 's/^"//' -e 's/,$//' -e 's/"$//' "$FILELIST" \
                       | grep -v '^[[:space:]]*$' | grep -v '^#')

NFILES=${#INFILES[@]}
[[ $NFILES -gt 0 ]] || { echo "ERROR: la lista '$FILELIST' esta vacia."; exit 1; }

echo "=============================================="
echo " Macro      : $MACRO_ABS"
echo " Lista      : $FILELIST  ($NFILES ficheros)"
echo " Salidas    : $OUTDIR"
echo " Logs       : $LOGDIR"
echo " Fusionado  : $MERGED"
echo " Nucleos    : $NCORES"
echo "=============================================="

# ------------------------------- worker ----------------------------------
process_one()
{
    local infile="$1"
    local base out log rc

    base=$(basename "$infile" .root)
    out="${OUTDIR}/${base}_reduced.root"
    log="${LOGDIR}/${base}.log"

    if [[ ! -f "$infile" ]]; then
        echo "[FALTA ] $infile (no existe)"
        return 1
    fi

    if [[ -s "$out" && -z "${FORCE:-}" ]]; then
        echo "[SALTO ] $base (ya existe, usa FORCE=1 para rehacerlo)"
        return 0
    fi

    root -l -b -q "${MACRO_ABS}(\"${infile}\",\"${out}\")" > "$log" 2>&1
    rc=$?

    if [[ $rc -ne 0 || ! -s "$out" ]]; then
        echo "[FALLO ] $base  (rc=$rc)  -> $log"
        return 1
    fi

    echo "[OK    ] $base -> $out"
    return 0
}
export -f process_one
export OUTDIR LOGDIR MACRO_ABS FORCE

# ---------------------------- procesado paralelo -------------------------
if [[ -z "${MERGE_ONLY:-}" ]]; then
    START=$(date +%s)

    printf '%s\0' "${INFILES[@]}" \
        | xargs -0 -n1 -P "$NCORES" -I{} bash -c 'process_one "$@"' _ {}
    XRC=$?

    END=$(date +%s)
    echo "----------------------------------------------"
    echo "Procesado terminado en $((END-START)) s."
    if [[ $XRC -ne 0 ]]; then
        echo "AVISO: al menos un job fallo. Revisa $LOGDIR (grep -l ERROR $LOGDIR/*.log)."
    fi
fi

# --------------------------------- hadd ----------------------------------
if [[ -n "${NO_MERGE:-}" ]]; then
    echo "NO_MERGE=1 -> no se fusiona."
    exit 0
fi

shopt -s nullglob
OUTS=( "${OUTDIR}"/*_reduced.root )
shopt -u nullglob

if [[ ${#OUTS[@]} -eq 0 ]]; then
    echo "ERROR: no hay ficheros de salida que fusionar en $OUTDIR."
    exit 1
fi

echo "----------------------------------------------"
echo "Fusionando ${#OUTS[@]} ficheros en $MERGED ..."

# -f  sobrescribe la salida
# -k  ignora ficheros corruptos en vez de abortar
# -j  hadd en paralelo (ROOT >= 6.20); si no existe, se reintenta sin -j
if ! hadd -f -k -j "$NCORES" "$MERGED" "${OUTS[@]}"; then
    echo "hadd -j no disponible o fallo; reintentando en serie..."
    hadd -f -k "$MERGED" "${OUTS[@]}" || { echo "ERROR: hadd fallo."; exit 1; }
fi

echo "=============================================="
echo "Listo: $MERGED"
root -l -b -q -e "TFile::Open(\"${MERGED}\")->Get(\"outtree\")->Print(\"toponly\")" 2>/dev/null \
    || echo "(no se pudo imprimir el resumen del tree, pero el fichero esta creado)"