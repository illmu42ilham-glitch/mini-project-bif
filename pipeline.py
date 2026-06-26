import requests
import matplotlib
matplotlib.use('Agg')  # Mencegah GUI matplotlib muncul di server backend
import matplotlib.pyplot as plt
import io
import base64
import pandas as pd

# ==========================================
# POIN 1: INTEGRASI NCBI & PENYIMPANAN LIST
# ==========================================
def fetch_fasta_from_ncbi(search_term):
    """Mencari nama protein/organisme di NCBI dan mengambil teks FASTA-nya."""
    try:
        # Step 1: Cari ID Nukleotida berdasarkan kata kunci pengguna
        search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term={search_term}&retmode=json&retmax=5"
        search_res = requests.get(search_url).json()
        id_list = search_res.get("esearchresult", {}).get("idlist", [])
        
        if not id_list:
            return None
        
        # Step 2: Ambil data sekuens FASTA berdasarkan ID yang ditemukan [cite: 101]
        ids = ",".join(id_list)
        fetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={ids}&rettype=fasta&retmode=text"
        fasta_data = requests.get(fetch_url).text
        return fasta_data
    except Exception as e:
        print(f"Error fetching NCBI: {e}")
        return None

def parse_fasta(fasta_text):
    """Membaca teks FASTA dan menyimpannya ke dalam struktur data List[cite: 101, 102]."""
    sequences_list = []  # List untuk menampung objek sekuens [cite: 102]
    current_header = ""
    current_seq = []
    
    for line in fasta_text.strip().split("\n"):
        if line.startswith(">"):
            if current_header:
                sequences_list.append({"header": current_header, "sequence": "".join(current_seq)})
            current_header = line[1:]
            current_seq = []
        else:
            current_seq.append(line.strip().upper())
            
    if current_header:
        sequences_list.append({"header": current_header, "sequence": "".join(current_seq)})
        
    return sequences_list

# ==========================================
# POIN 2: ANALISIS DICTIONARY & GC CONTENT
# ==========================================
def analyze_sequences(sequences):
    """Menghitung frekuensi nukleotida dan persentase GC Content[cite: 103, 104]."""
    analyzed_data = []
    
    for item in sequences:
        header = item["header"]
        seq = item["sequence"]
        
        # Menggunakan Dictionary untuk menghitung frekuensi basa [cite: 103]
        freq = {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0}
        for base in seq:
            if base in freq:
                freq[base] += 1
            else:
                freq["N"] += 1  # Untuk karakter tak dikenal
                
        total_bases = sum(freq.values())
        
        # Menghitung nilai GC Content
        gc_content = ((freq["G"] + freq["C"]) / total_bases * 100) if total_bases > 0 else 0
        
        analyzed_data.append({
            "short_header": header[:40] + "...",  # Potongan teks untuk visualisasi halaman web
            "full_header": header,
            "frequencies": freq,
            "gc_content": round(gc_content, 2),
            "length": total_bases
        })
    
    # Mengurutkan sekuens berdasarkan GC Content tertinggi secara descending [cite: 104]
    analyzed_data.sort(key=lambda x: x["gc_content"], reverse=True)
    return analyzed_data

# ==========================================
# POIN 3: VISUALISASI GRAFIK HASIL GC
# ==========================================
def generate_gc_chart(analyzed_data):
    """Membuat grafik batang GC Content menggunakan Matplotlib dan diubah ke base64[cite: 106]."""
    headers = [d["short_header"][:12] for d in analyzed_data]
    gc_values = [d["gc_content"] for d in analyzed_data]
    
    plt.figure(figsize=(7, 4))
    bars = plt.bar(headers, gc_values, color='#0ea5e9', edgecolor='#0284c7', linewidth=1.2)
    
    # Mempercantik tampilan grafik
    plt.ylabel('GC Content (%)', fontsize=11, fontweight='bold')
    plt.title('Grafik Analisis GC Content', fontsize=13, fontweight='bold', pad=15)
    plt.ylim(0, 100)
    plt.xticks(rotation=15, ha='right')
    
    # Menambahkan label nilai di atas batang grafik
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2.0, yval + 2, f'{yval}%', ha='center', va='bottom', fontsize=9)
        
    plt.tight_layout()
    
    # Simpan grafik ke memori agar bisa langsung dibaca oleh HTML tanpa membuat file sampah
    img = io.BytesIO()
    plt.savefig(img, format='png', bbox_inches='tight')
    img.seek(0)
    plot_url = base64.b64encode(img.getvalue()).decode('utf8')
    plt.close()
    return plot_url

# ==========================================
# POIN 4: MENULIS DATA KE CSV
# ==========================================
def export_to_csv(analyzed_data, filepath):
    """Menuliskan hasil akhir pipeline ke dalam file CSV menggunakan Pandas[cite: 107]."""
    rows = []
    for d in analyzed_data:
        rows.append({
            "Header": d["full_header"],
            "Length": d["length"],
            "GC_Content(%)": d["gc_content"],
            "A": d["frequencies"]["A"],
            "T": d["frequencies"]["T"],
            "G": d["frequencies"]["G"],
            "C": d["frequencies"]["C"],
            "N": d["frequencies"]["N"]
        })
    df = pd.DataFrame(rows)
    df.to_csv(filepath, index=False)