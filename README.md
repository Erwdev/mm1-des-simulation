# mm1-des-simulation

# Discrete Event Simulation (DES) M/M/1 Queue System

Simulasi sistem antrian M/M/1 menggunakan C++ dengan dukungan multi-replikasi, confidence interval, dan berbagai mode terminasi.

## 📋 Deskripsi Proyek

Program ini mengimplementasikan Discrete Event Simulation untuk model antrian M/M/1 (Markovian arrival, Markovian service, 1 server) dengan fitur-fitur:

- **Dua Mode Terminasi**: Berdasarkan jumlah pelanggan terlayani atau waktu simulasi
- **Multi-Replikasi**: Mendukung multiple runs untuk analisis statistik
- **Confidence Interval**: Perhitungan CI 95% otomatis
- **Warm-up Period**: Eliminasi bias transient state
- **Flexible Scenarios**: Dukungan untuk variasi eksperimen (λ(t), finite queue, deterministic service)
- **CSV Export**: Output terstruktur untuk analisis lanjutan

## 🎓 Anggota Kelompok

- [Nama Anggota 1] - [NIM] - [Role: Coding/Testing]
- [Nama Anggota 2] - [NIM] - [Role: Analisis Data]
- [Nama Anggota 3] - [NIM] - [Role: Dokumentasi]
- [Nama Anggota 4] - [NIM] - [Role: Presentasi]

## 🛠️ Requirements

### Software Requirements
- **Compiler**: g++ atau clang++ dengan support C++17 atau lebih baru
- **Operating System**: Linux, macOS, atau Windows (dengan MinGW/WSL)
- **Memory**: Minimal 512 MB RAM
- **Storage**: ~10 MB untuk program dan output files

### Library Dependencies
Program ini menggunakan **standard library C++ only**, tidak memerlukan library eksternal:
- `<queue>` - untuk priority queue (FEL)
- `<random>` - untuk random number generation
- `<fstream>` - untuk CSV output
- `<cmath>` - untuk perhitungan statistik
- `<vector>` - untuk storage hasil replikasi

## 📦 Struktur Proyek

```
des_mm1_simulator/
├── README.md                    # File ini
├── des_sim.cpp                  # Source code utama
├── Makefile                     # Build automation (optional)
├── docs/
│   ├── laporan.pdf              # Laporan praktikum
│   └── presentasi.pdf           # Slide presentasi
├── output/
│   ├── results_per_rep.csv      # Hasil per replikasi
│   └── summary.csv              # Ringkasan statistik
└── tests/
    └── test_scenarios.sh        # Script untuk testing otomatis
```

## 🔧 Instalasi & Kompilasi

### Kompilasi dengan g++

```bash
# Kompilasi program
g++ -std=c++17 -O2 -Wall -o des_sim des_sim.cpp

# Atau dengan optimasi maksimal
g++ -std=c++17 -O3 -march=native -Wall -o des_sim des_sim.cpp
```

### Kompilasi dengan clang++

```bash
clang++ -std=c++17 -O2 -Wall -o des_sim des_sim.cpp
```

### Menggunakan Makefile (jika tersedia)

```bash
make          # Kompilasi program
make clean    # Hapus binary dan output files
make test     # Jalankan test scenarios
```

## 🚀 Cara Penggunaan

### Sintaks Umum

```bash
./des_sim [OPTIONS]
```

### Parameter Wajib

| Parameter | Deskripsi | Tipe | Default |
|-----------|-----------|------|---------|
| `--lambda` | Arrival rate (λ) | float | 0.9 |
| `--mu` | Service rate (μ) | float | 1.0 |
| `--term` | Mode terminasi: `served` atau `time` | string | served |
| `--reps` | Jumlah replikasi | int | 10 |

### Parameter Opsional

| Parameter | Deskripsi | Tipe | Default |
|-----------|-----------|------|---------|
| `--maxServed` | Jumlah pelanggan (mode `served`) | int | 10000 |
| `--horizonT` | Waktu simulasi (mode `time`) | float | 10000.0 |
| `--warmup` | Periode warm-up | int | 1000 |
| `--seed` | Random seed | int | 12345 |
| `--queueCap` | Kapasitas antrian (-1 = unlimited) | int | -1 |
| `--outdir` | Direktori output | string | ./ |
| `--help` | Tampilkan bantuan | - | - |

## 📝 Contoh Penggunaan

### 1. Simulasi Standar M/M/1 (Terminasi by Served)

```bash
./des_sim --lambda 0.9 --mu 1.0 --maxServed 20000 --warmup 1000 --reps 10 --seed 12345 --term served --outdir ./output/
```

**Deskripsi**: Simulasi dengan ρ = 0.9 (heavy traffic), terminasi setelah 20,000 pelanggan terlayani, 10 replikasi.

### 2. Simulasi dengan Terminasi by Time

```bash
./des_sim --lambda 0.5 --mu 1.0 --horizonT 50000 --warmup 1000 --reps 20 --seed 99999 --term time --outdir ./output/
```

**Deskripsi**: Simulasi dengan ρ = 0.5 (moderate traffic), terminasi setelah t = 50,000, 20 replikasi.

### 3. Simulasi Light Traffic

```bash
./des_sim --lambda 0.1 --mu 1.0 --maxServed 10000 --warmup 500 --reps 15 --term served
```

**Deskripsi**: Simulasi dengan ρ = 0.1 (light traffic), cocok untuk validasi teori.

### 4. Simulasi dengan Finite Queue

```bash
./des_sim --lambda 0.9 --mu 1.0 --maxServed 15000 --queueCap 50 --warmup 1000 --reps 10 --term served
```

**Deskripsi**: Simulasi dengan kapasitas antrian terbatas (maksimal 50 pelanggan dalam sistem).

### 5. Batch Testing Multiple Scenarios

```bash
# Light traffic
./des_sim --lambda 0.3 --mu 1.0 --maxServed 10000 --reps 10 --term served --outdir ./output/light/

# Moderate traffic
./des_sim --lambda 0.5 --mu 1.0 --maxServed 10000 --reps 10 --term served --outdir ./output/moderate/

# Heavy traffic
./des_sim --lambda 0.9 --mu 1.0 --maxServed 10000 --reps 10 --term served --outdir ./output/heavy/
```

## 📊 Format Output

### 1. results_per_rep.csv

File ini berisi hasil dari setiap replikasi individual.

```csv
RepID,AvgQ,Utilization,AvgDelay,AvgWait,NumServed,SimTime
1,8.1234,0.8998,9.0123,8.0345,20000,22234.56
2,8.3456,0.9012,9.2345,8.2567,20000,22145.78
...
```

**Kolom**:
- `RepID`: Nomor replikasi (1, 2, 3, ...)
- `AvgQ`: Rata-rata jumlah pelanggan dalam sistem (L)
- `Utilization`: Utilisasi server (ρ)
- `AvgDelay`: Rata-rata waktu di sistem (W)
- `AvgWait`: Rata-rata waktu menunggu di antrian (Wq)
- `NumServed`: Jumlah pelanggan yang terlayani
- `SimTime`: Total waktu simulasi

### 2. summary.csv

File ini berisi ringkasan statistik dengan confidence interval.

```csv
Metric,Mean,StdDev,CI_Lower,CI_Upper,CI_Width
AvgQ,8.2345,0.1234,8.1456,8.3234,0.1778
Utilization,0.9005,0.0012,0.8998,0.9012,0.0014
AvgDelay,9.1234,0.1456,9.0123,9.2345,0.2222
AvgWait,8.1345,0.1356,8.0234,8.2456,0.2222
```

**Kolom**:
- `Metric`: Nama metrik performa
- `Mean`: Nilai rata-rata dari semua replikasi
- `StdDev`: Standard deviasi
- `CI_Lower`: Batas bawah confidence interval 95%
- `CI_Upper`: Batas atas confidence interval 95%
- `CI_Width`: Lebar interval kepercayaan

## 📈 Interpretasi Hasil

### Metrik Utama

1. **Average Queue Length (L atau AvgQ)**
   - Rata-rata jumlah pelanggan dalam sistem (menunggu + dilayani)
   - Teoritis M/M/1: L = ρ / (1 - ρ)
   - Semakin tinggi ρ, semakin panjang antrian

2. **Utilization (ρ)**
   - Proporsi waktu server sibuk
   - ρ = λ / μ
   - Sistem stabil jika ρ < 1

3. **Average Delay (W atau AvgDelay)**
   - Rata-rata waktu pelanggan di sistem (waiting + service)
   - Teoritis M/M/1: W = 1 / (μ - λ)
   - Meningkat drastis saat ρ mendekati 1

4. **Average Wait (Wq atau AvgWait)**
   - Rata-rata waktu menunggu di antrian (tidak termasuk service)
   - Teoritis M/M/1: Wq = ρ / (μ - λ)

### Validasi dengan Little's Law

Verifikasi hasil simulasi menggunakan Little's Law:
- **L = λ × W** (Jumlah dalam sistem = arrival rate × waktu di sistem)
- **Lq = λ × Wq** (Jumlah dalam antrian = arrival rate × waktu menunggu)

### Analisis Confidence Interval

- **CI Width sempit** → Hasil stabil, replikasi cukup
- **CI Width lebar** → Perlu lebih banyak replikasi
- **CI tidak overlap** → Perbedaan signifikan antar skenario

## 🧪 Validasi & Testing

### Test Case 1: Light Traffic (ρ = 0.3)

```bash
./des_sim --lambda 0.3 --mu 1.0 --maxServed 10000 --warmup 500 --reps 10 --term served
```

**Expected Results**:
- L ≈ 0.43
- W ≈ 1.43
- Utilization ≈ 0.30

### Test Case 2: Moderate Traffic (ρ = 0.5)

```bash
./des_sim --lambda 0.5 --mu 1.0 --maxServed 10000 --warmup 1000 --reps 10 --term served
```

**Expected Results**:
- L ≈ 1.0
- W ≈ 2.0
- Utilization ≈ 0.50

### Test Case 3: Heavy Traffic (ρ = 0.9)

```bash
./des_sim --lambda 0.9 --mu 1.0 --maxServed 20000 --warmup 2000 --reps 15 --term served
```

**Expected Results**:
- L ≈ 9.0
- W ≈ 10.0
- Utilization ≈ 0.90

### Automated Testing

```bash
# Jika test script tersedia
chmod +x tests/test_scenarios.sh
./tests/test_scenarios.sh
```

## 🔬 Skenario Eksperimen Lanjutan

### Scenario 1: Time-Varying Arrival Rate λ(t)

Implementasi arrival rate yang berubah berdasarkan waktu (misal: pagi sibuk, siang sepi).

```cpp
// Dalam kode: modifikasi fungsi lambda
double getLambda(double t) {
    // Contoh: sinusoidal pattern
    return baseλ + amplitudo * sin(2 * π * t / period);
}
```

### Scenario 2: Finite Queue Capacity

Sistem dengan kapasitas antrian terbatas (M/M/1/K).

```bash
./des_sim --lambda 0.9 --mu 1.0 --maxServed 15000 --queueCap 10 --reps 10 --term served
```

**Analisis**: Hitung probability of blocking = (jumlah ditolak) / (jumlah arrival)

### Scenario 3: Deterministic Service Time (M/D/1)

Service time konstan D = 1/μ (bukan eksponensial).

```cpp
// Dalam kode: ganti exponential dengan konstan
double serviceTime = 1.0 / mu;  // deterministic
```

**Expected**: Performa lebih baik dari M/M/1 (variance lebih rendah)

## ⚠️ Troubleshooting

### Error: "Command not found: ./des_sim"

**Solusi**:
```bash
# Pastikan file executable
chmod +x des_sim

# Atau compile ulang
g++ -std=c++17 -O2 -o des_sim des_sim.cpp
```

### Error: Compilation failed

**Solusi**:
```bash
# Check compiler version (harus support C++17)
g++ --version

# Coba compile dengan verbose
g++ -std=c++17 -O2 -Wall -v -o des_sim des_sim.cpp
```

### Warning: "System is unstable (λ ≥ μ)"

**Penyebab**: ρ ≥ 1.0, sistem tidak stabil

**Solusi**: Pastikan λ < μ untuk sistem stabil

### Output: Hasil tidak match teori

**Possible Issues**:
1. Warm-up period terlalu pendek → Tingkatkan `--warmup`
2. Jumlah replikasi terlalu sedikit → Tingkatkan `--reps`
3. Bug dalam implementasi → Review logika DES

### Performance: Simulasi terlalu lambat

**Solusi**:
```bash
# Compile dengan optimasi maksimal
g++ -std=c++17 -O3 -march=native -o des_sim des_sim.cpp

# Kurangi jumlah pelanggan atau waktu simulasi
./des_sim --lambda 0.9 --mu 1.0 --maxServed 5000 --reps 5
```

## 📚 Referensi

1. **Law, A. M.** (2015). *Simulation Modeling and Analysis* (5th ed.). McGraw-Hill.
2. **Ross, S. M.** (2014). *Introduction to Probability Models* (11th ed.). Academic Press.
3. **Banks, J., et al.** (2010). *Discrete-Event System Simulation* (5th ed.). Pearson.
4. **Gross, D., & Harris, C. M.** (1998). *Fundamentals of Queueing Theory* (3rd ed.). Wiley.

### Online Resources

- [Queueing Theory Overview](https://en.wikipedia.org/wiki/Queueing_theory)
- [M/M/1 Queue Model](https://en.wikipedia.org/wiki/M/M/1_queue)
- [Little's Law Explanation](https://en.wikipedia.org/wiki/Little%27s_law)
- [C++ Priority Queue Documentation](https://en.cppreference.com/w/cpp/container/priority_queue)

## 📧 Kontak

Untuk pertanyaan atau issue terkait program ini:

- **Email Kelompok**: [email@example.com]
- **GitHub Repository**: [https://github.com/username/des-mm1-simulator]
- **Dosen Pengampu**: [Nama Dosen] - [email dosen]

## 📄 Lisensi

Program ini dibuat untuk keperluan akademis Praktikum [Nama Mata Kuliah] - [Semester/Tahun].

© 2025 [Nama Kelompok]. All rights reserved.

---

**Last Updated**: [Tanggal]  
**Version**: 1.0.0  
**Status**: ✅ Production Ready