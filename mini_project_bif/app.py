from flask import Flask, render_template, request, send_file
import os
from pipeline import fetch_fasta_from_ncbi, parse_fasta, analyze_sequences, generate_gc_chart, export_to_csv

app = Flask(__name__)
OUTPUT_DIR = "outputs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

@app.route("/", methods=["GET", "POST"])
def index():
    results = None
    top_3 = None
    chart_url = None
    error = None
    search_query = ""

    if request.method == "POST":
        search_query = request.form.get("search_query")
        if search_query:
            # Jalankan integrasi pencarian NCBI API [cite: 101]
            fasta_text = fetch_fasta_from_ncbi(search_query)
            
            if fasta_text:
                # Konversi teks ke struktur List [cite: 102]
                raw_sequences = parse_fasta(fasta_text)
                
                if raw_sequences:
                    # Proses analisis Dictionary & GC sorting [cite: 103, 104]
                    results = analyze_sequences(raw_sequences)
                    
                    # Ambil 3 sekuens terbaik [cite: 105]
                    top_3 = results[:3]
                    
                    # Buat grafik visualisasi [cite: 106]
                    chart_url = generate_gc_chart(results)
                    
                    # Tulis hasil ke CSV secara otomatis [cite: 107]
                    csv_path = os.path.join(OUTPUT_DIR, "hasil_analisis.csv")
                    export_to_csv(results, csv_path)
                else:
                    error = "Gagal memproses parsing data FASTA dari sekuens tersebut."
            else:
                error = f"Protein atau organisme '{search_query}' tidak ditemukan di database NCBI Nuccore."
                
    return render_template("index.html", results=results, top_3=top_3, chart_url=chart_url, error=error, query=search_query)

@app.route("/download")
def download():
    """Route untuk mengunduh berkas CSV hasil ekspor[cite: 107]."""
    csv_path = os.path.join(OUTPUT_DIR, "hasil_analisis.csv")
    if os.path.exists(csv_path):
        return send_file(csv_path, as_attachment=True)
    return "File tidak ditemukan.", 404

if __name__ == "__main__":
    app.run(debug=True)