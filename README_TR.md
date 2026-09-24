# AI Plant Gene Analyzer

## Uygulamanın Mevcut İşlevleri

Bu belge, uygulamanın 24 Eylül 2026 tarihindeki mevcut işlevsel durumunu özetler. Uygulama, DNA, RNA ve protein dizilerinin analizinden başlayarak benzerlik araması, varyant analizi, çeviri, filogeni ve biyokimyasal incelemeye kadar birçok işlemi tek bir Streamlit arayüzünde sunar.

## Uygulamayı Başlatma

Proje klasöründe aşağıdaki komut çalıştırılabilir:

```powershell
.\.venv\Scripts\python.exe -m streamlit run app.py
```

Uygulama iki ana sayfadan oluşur:

- **Sequence Analysis**: Bir veya birden fazla dizinin kapsamlı analizi.
- **Independent Tools**: Ana analiz akışından bağımsız olarak kullanılabilen araçlar.

# 1. Sequence Analysis

## Girdi Türleri

Kullanıcı aşağıdaki yöntemlerle veri sağlayabilir:

- DNA dizisi yapıştırma
- RNA dizisi yapıştırma
- Protein dizisi yapıştırma
- `.fasta`, `.fa` veya `.txt` dosyası yükleme
- Aynı FASTA dosyasında birden fazla dizi kullanma
- İsteğe bağlı referans dizisi ekleme
- Hazır demo dizilerinden birini seçme

Dizi tipi otomatik olarak algılanabilir veya kullanıcı tarafından DNA, RNA ya da protein olarak seçilebilir.

Birden fazla FASTA kaydı yüklendiğinde:

- Tek bir kayıt seçilebilir.
- Tüm kayıtlar toplu olarak analiz edilebilir.
- Toplu analiz sonunda her dizi için uzunluk, tip, en iyi eşleşme ve benzerlik özeti gösterilir.

## İstatistiksel Analiz

### DNA ve RNA dizileri için

- Dizi uzunluğu
- GC ve AT yüzdesi
- GC/AT oranı
- Kayan pencere ile GC profili
- GC skew profili
- A, T, G ve C sayıları
- Başlangıç ve dur kodonu kontrolü
- Tam ORF kontrolü
- Kodon kullanımı
- Codon Adaptation Index (CAI), uygun organizma referansı varsa
- CG, CHG ve CHH metilasyon bağlamları
- Düşük karmaşıklık veya tekrarlı bölgeler
- Restriksiyon enzimi bölgeleri
- Primer tasarımı için başlangıç önerileri
- Bilinen düzenleyici motifler
- Organizma referansına göre uzunluk ve GC karşılaştırması

### Protein dizileri için

- Amino asit uzunluğu
- Amino asidi bileşimi
- Benzersiz rezidü sayısı
- En sık bulunan rezidü
- Moleküler ağırlık
- İzolektrik nokta
- Hidrofobiklik
- GRAVY değeri
- Kararsızlık indeksi
- Alifatik indeks
- Bilinen protein motifleri

Sonuçlar tablolar, metrikler ve etkileşimli grafiklerle gösterilir.

## Benzerlik Araması

Girilen dizi, yerel gen veritabanındaki kayıtlarla karşılaştırılır.

Gösterilen bilgiler:

- En iyi eşleşen genler
- Benzerlik yüzdesi
- Global hizalama skoru
- Eşleşmeler ve uyumsuzluklar
- Gap sayısı
- Global ve lokal kapsama
- Organizmalar
- Gen sembolleri ve accession bilgileri
- Özellik veya trait bilgileri
- Gen açıklamaları
- Veri kaynağı
- Benzerlik güven sınıfı

İki arama modu vardır:

- **Balanced Search**: Daha hızlıdır; adayları önceden filtreler.
- **Deep Search**: Daha geniş ve hassas arama yapar, ancak daha uzun sürer.

Ek görseller ve karşılaştırmalar:

- Benzerlik skor grafiği
- En iyi üç adayın karşılaştırma tablosu
- Hizalama görünümü
- Kapsama heatmap'i
- Güven göstergesi
- En iyi adaylar arasında çoklu hizalama
- Korunmuş ve değişken bölgelerin görünümü

Benzerlik sonucu, tek başına gen kimliği veya biyolojik fonksiyon kanıtı olarak kabul edilmemelidir.

## Mutasyon ve Varyant Analizi

Kullanıcı bir referans dizi sağlarsa, girilen dizi bu referansla karşılaştırılır. Açık bir referans verilmezse, yeterince yüksek benzerlik bulunduğunda en iyi veritabanı eşleşmesi referans olarak kullanılabilir.

Tespit edilen bilgiler:

- Substitüsyonlar
- Transition ve transversion olayları
- Inserisyonlar ve delesyonlar
- Mutasyon oranı
- Gap'li ve gap'siz kimlik
- Variant konumları
- Olası sonuçlar
- Sessiz mutasyonlar
- Missense ve nonsense değişiklikler
- Readthrough olayları
- Frameshift olayları
- Potansiyel olarak önemli varyantlar

Sunulan görünümler:

- Mutasyon haritası
- Pencere bazında varyant frekansı
- Dizi-referans hizalaması
- Sadece önemli mutasyonları gösterme filtresi

Dışa aktarma biçimleri:

- Mutasyon CSV dosyası
- Substitüsyon VCF dosyası

## Çeviri ve ORF Analizi

DNA veya RNA dizileri için uygulama:

- `+1`, `+2`, `+3`, `-1`, `-2`, `-3` okuma çerçevelerini analiz eder.
- Seçilen çerçeveyi proteine çevirir.
- Dur kodonunu ve tam ORF durumunu kontrol eder.
- Protein uzunluğunu hesaplar.
- Belirsiz nükleotidleri işaretler.
- Kodon haritası oluşturur.
- Altı okuma çerçevesini karşılaştırır.
- İnceleme için önerilen çerçeveyi belirtir.
- ORF'leri bulur ve grafik üzerinde gösterir.
- Tamamlayıcı ve ters tamamlayıcı dizileri üretir.

Dışa aktarma biçimleri:

- Çevrilmiş protein FASTA dosyası
- Tahmin edilen ORF'ler için GFF3 dosyası

Protein girdilerinde çeviri yapılmaz; bunun yerine protein özellikleri gösterilir.

## Yapay Zekâ Destekli Yorumlama

Uygulama analiz sonuçlarını özetleyen bir yorum üretir. Bu yorumda aşağıdaki noktalar ele alınabilir:

- Dizi tipi ve uzunluğu
- GC oranı veya protein özellikleri
- ORF durumu
- Dizi kalite bilgileri
- Motifler
- Benzerlik seviyesi
- İncelenmesi gereken noktalar

Bu yorum biyolojik doğrulamanın veya deneysel doğrulamanın yerine geçmez.

## Ham Dizi Görünümü

Uygulama tarafından temizlenerek analiz edilen dizi ayrıca görüntülenebilir. Böylece:

- Kullanılan gerçek dizi kontrol edilebilir.
- Boşlukların ve gereksiz karakterlerin kaldırıldığı görülebilir.
- Analiz aşamalarına gönderilen veri doğrulanabilir.

## Sonuçları Dışa Aktarma

Ana analiz sayfasında kullanılabilen rapor biçimleri arasında şunlar bulunur:

- JSON
- CSV
- HTML
- Mutasyon CSV
- VCF
- FASTA
- GFF3

# 2. Independent Tools

Bu sayfa, ana dizi analizini çalıştırmadan belirli araçların ayrı ayrı kullanılmasını sağlar.

## Hizalamalar

### Çoklu dizi hizalaması (MSA)

İki veya daha fazla homolog dizi kullanılabilir.

Özellikler:

- DNA veya protein dizileri
- Otomatik dizi tipi algılama
- Referans dizi seçimi
- Protein substitüsyon matrisi seçimi
- Gap açma ve uzatma cezaları
- Hizalanmış kolon sayısı
- Korunma yüzdesi
- Konsensus dizisi
- Değişken kolonlar
- Düşük karmaşıklık uyarısı

Dışa aktarma biçimleri:

- Aligned FASTA
- Clustal ALN
- Metrics CSV
- PHYLIP
- NEXUS

### İkili hizalama

Tam olarak iki dizi için:

- Global Needleman-Wunsch hizalaması
- Lokal Smith-Waterman hizalaması
- Kimlik yüzdesi
- Eşleşme sayısı
- Uyumsuzluk sayısı
- Gap sayısı
- Hizalama skoru
- Her iki dizi için kapsama
- Gap'ler hariç kimlik

DNA-DNA veya protein-protein karşılaştırması yapılabilir.

## Uzaklık Matrisi

Birden fazla dizi kullanılarak uzaklık matrisi hesaplanabilir.

Mevcut yöntemler:

- Hamming
- Jukes-Cantor
- Kimura
- Proteinler için PAM

Sonuçlar:

- Uzaklık heatmap'i
- Sayısal uzaklık matrisi
- En yakın dizi çifti
- En uzak dizi çifti
- CSV indirme seçeneği

Hesaplanan matris doğrudan filogenetik ağaç oluşturma aracına aktarılabilir.

## Filogenetik Ağaç

Ağaç oluşturma yöntemleri:

- UPGMA
- Neighbor Joining

Gösterilen sonuçlar:

- UPGMA için dendrogram
- Neighbor Joining için ağaç grafiği
- Kullanılan algoritmanın bilgisi
- Kullanılan uzaklık yöntemi
- Newick formatında dışa aktarma

Ağaç, hesaplanan dizisel uzaklıklara dayalı bir gruplandırmadır; tek başına biyolojik soy ilişkisini kanıtlamaz.

## Protein Biyokimyasal Analizi

Bir protein dizisi girildiğinde:

- Uzunluk
- Moleküler ağırlık
- İzolektrik nokta
- Hidrofobiklik
- Amino asit dağılımı

hesaplanır ve amino asit bileşimi grafikle gösterilir.

## Trait Araması

Toplanan gen ve tür verileri arasında trait/özellik araması yapılabilir. Bu bölüm, doğrudan dizi analizinden bağımsız olarak genlerle ve türlerle ilişkili biyolojik bilgilerin incelenmesini sağlar.

# Genel Değerlendirme

Uygulamanın mevcut sürümü şu işlemleri tek bir arayüzde birleştirir:

- Dizi temizleme ve doğrulama
- DNA, RNA ve protein istatistikleri
- GC, kodon ve metilasyon analizi
- Protein özellikleri
- Benzerlik araması
- Global, lokal ve çoklu hizalama
- Mutasyon ve varyant analizi
- Çeviri ve ORF tahmini
- Uzaklık matrisi
- Filogenetik ağaç oluşturma
- Trait araması
- Etkileşimli grafikler ve tablolar
- JSON, CSV, HTML, FASTA, VCF, GFF3, PHYLIP, NEXUS, Clustal ve Newick dışa aktarımları

Bu nedenle uygulama, ham bir diziden başlayarak istatistiksel analiz, gen benzerliği, varyant incelemesi, anotasyon ve filogenetik değerlendirmeye kadar uzanan bütünleşik bir biyoinformatik analiz platformu olarak kullanılabilir.
