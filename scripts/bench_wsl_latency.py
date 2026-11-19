"""Mesure ou le temps est perdu (WSL, DIAMOND proteome, DIAMOND Swiss-Prot, hmmscan). Lecture seule.
Usage : python scripts\\bench_wsl_latency.py results\\NP_001311551.2.fasta"""
import subprocess, sys, time
from pathlib import Path

REMOTE = r'''
R="$HOME/reference"; W=$(mktemp -d /tmp/bench.XXXXXX); cat > "$W/q.faa"
TIMEFORMAT=%R
echo "== MACHINE =="; echo "nproc=$(nproc)"; free -m | sed -n '1,3p'
echo "reference sur: $(df -P "$R" | tail -1 | awk '{print $1, $6}')"
ls -la "$R" | awk '{print $5, $9}' | sed -n '1,40p'
cpu() { echo "$1 : $( { time "${@:2}" >/dev/null 2>&1; } 2>&1 ) s"; }
N=$(nproc)
for T in 2 $N; do
  echo "== threads=$T =="
  cpu "proteome sensitive (1er appel)" diamond blastp -d "$R/plant_proteomes" -q "$W/q.faa" --sensitive --evalue 1e-5 -k 45 --threads $T -o "$W/o1.tsv"
  cpu "proteome sensitive (2e appel)"  diamond blastp -d "$R/plant_proteomes" -q "$W/q.faa" --sensitive --evalue 1e-5 -k 45 --threads $T -o "$W/o1.tsv"
  cpu "proteome fast"                  diamond blastp -d "$R/plant_proteomes" -q "$W/q.faa" --fast --evalue 1e-5 -k 45 --threads $T -o "$W/o1.tsv"
  cpu "swissprot sensitive (1er)"      diamond blastp -d "$R/swissprot_diamond" -q "$W/q.faa" --sensitive --evalue 1e-5 --max-target-seqs 5 --threads $T -o "$W/o2.tsv"
  cpu "swissprot sensitive (2e)"       diamond blastp -d "$R/swissprot_diamond" -q "$W/q.faa" --sensitive --evalue 1e-5 --max-target-seqs 5 --threads $T -o "$W/o2.tsv"
  cpu "hmmscan Pfam (1er)"             hmmscan --cut_ga --noali --cpu $T --domtblout "$W/d.tbl" "$R/Pfam-A.hmm" "$W/q.faa"
  cpu "hmmscan Pfam (2e)"              hmmscan --cut_ga --noali --cpu $T --domtblout "$W/d.tbl" "$R/Pfam-A.hmm" "$W/q.faa"
done
rm -rf "$W"; echo FIN
'''

def main():
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    fasta = Path(sys.argv[1]).read_bytes().replace(b"\r\n", b"\n")
    t = []
    for _ in range(3):
        t0 = time.time()
        subprocess.run(["wsl.exe", "-d", "Ubuntu", "--exec", "true"], capture_output=True)
        t.append(time.time() - t0)
    print("Demarrage WSL (3 essais, s) :", ", ".join(f"{x:.2f}" for x in t))
    t0 = time.time()
    r = subprocess.run(["wsl.exe", "-d", "Ubuntu", "--exec", "bash", "-c", REMOTE.replace("\r", "")],
                       input=fasta, capture_output=True, timeout=3000)
    print(r.stdout.decode("utf-8", "replace"))
    if r.stderr:
        print("STDERR:", r.stderr.decode("utf-8", "replace")[:500])
    print(f"Duree totale du banc : {time.time()-t0:.0f} s")

main()
