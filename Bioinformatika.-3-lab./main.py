import matplotlib.pyplot as plt
from bioinfokit.analys import fastq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


def ketvirto_uzd_a():
    for formatas in formatai:
        simbolis_uzkodavime = True
        for simbolis in skirtingi_simboliai:
            if ord(simbolis) < formatai[formatas][0] or ord(simbolis) > formatai[formatas][1]:
                simbolis_uzkodavime = False
        if simbolis_uzkodavime:
            print(formatas)


formatai = dict({
    'Sanger Phred+33': (33, 73),
    'Solexa Solexa+64': (59, 104),
    'Illumina 1.3+ Phred+64': (64, 104),
    'Illumina 1.5+ Phred+64': (67, 105),
    'Illumina 1.8+ Phred+33': (33, 74),
})


def ketvirto_uzd_b():
    for i in range(0, 100):
        abscises.append(i / 100)
    for x in abscises:
        ordinates.append(cg_dalmenys.count(x))
    plt.plot(abscises, ordinates)
    plt.xlabel('C/G nukletidų dalis read’o sekoje')
    plt.ylabel('read’ų skaičius')
    plt.show()


def ketvirto_uzd_c():
    # pirmas pikas 0.34
    print(abscises[ordinates[:40].index(max(ordinates[:40]))])
    # antras pikas 0.54
    print(abscises[40:60][ordinates[40:60].index(max(ordinates[40:60]))])
    # trecias pikas 0.70
    print(abscises[60:][ordinates[60:].index(max(ordinates[60:]))])

    pirmo_piko_ids = []
    pirmo_piko_sekos = []
    antro_piko_ids = []
    antro_piko_sekos = []
    trecio_piko_ids = []
    trecio_piko_sekos = []
    pirmo_piko_pozicijos = [j for j, x in enumerate(cg_dalmenys) if x == 0.34][:5]
    antro_piko_pozicijos = [j for j, x in enumerate(cg_dalmenys) if x == 0.54][:5]
    trecio_piko_pozicijos = [j for j, x in enumerate(cg_dalmenys) if x == 0.7][:5]
    for position in pirmo_piko_pozicijos:
        pirmo_piko_ids.append(ids[position])
        pirmo_piko_sekos.append(sekos[position])
    for position in antro_piko_pozicijos:
        antro_piko_ids.append(ids[position])
        antro_piko_sekos.append(sekos[position])
    for position in trecio_piko_pozicijos:
        trecio_piko_ids.append(ids[position])
        trecio_piko_sekos.append(sekos[position])
    print(pirmo_piko_sekos)
    print(antro_piko_sekos)
    print(trecio_piko_sekos)

    ids_sulpelis = []
    bakterijos_stulpelis = []
    for sekos_id in pirmo_piko_ids:
        ids_sulpelis.append(sekos_id)
    for sekos_id in antro_piko_ids:
        ids_sulpelis.append(sekos_id)
    for sekos_id in trecio_piko_ids:
        ids_sulpelis.append(sekos_id)

    for s in pirmo_piko_sekos:
        bakterijos_stulpelis.append(sprogimo_paieska(s))
    for s in antro_piko_sekos:
        bakterijos_stulpelis.append(sprogimo_paieska(s))
    for s in trecio_piko_sekos:
        bakterijos_stulpelis.append(sprogimo_paieska(s))

    print(ids_sulpelis)
    print(bakterijos_stulpelis)


def sprogimo_paieska(seq):
    rezultato_apdorojimas = NCBIWWW.qblast("blastn", "nt", seq, alignments=1, hitlist_size=1)
    irasai = NCBIXML.parse(rezultato_apdorojimas)
    for irasas in irasai:
        for lygiavimas in irasas.alignments:
            return lygiavimas.hit_def


if __name__ == '__main__':
    fastq_kartojimas = fastq.fastq_reader(file='duomenys\\reads_for_analysis.fastq')

    kokybes_ivertinimu_sujungimas = ""
    cg_dalmenys = []

    ids = []
    sekos = []

    for irasas in fastq_kartojimas:
        sekos_num, seka, _, kokybes_ivertinimas = irasas
        kokybes_ivertinimu_sujungimas = kokybes_ivertinimu_sujungimas + kokybes_ivertinimas
        cg_dalmenys.append(round((seka.count('G') + seka.count('C')) / len(seka), 2))
        ids.append(sekos_num)
        sekos.append(seka)
    skirtingi_simboliai = set(kokybes_ivertinimu_sujungimas)

    ketvirto_uzd_a()

    abscises = []
    ordinates = []
    ketvirto_uzd_b()

    ketvirto_uzd_c()
