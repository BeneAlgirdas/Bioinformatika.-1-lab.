# from Bio.Seq import Seq
# my_seq = Seq("AGTACACTGGT")
# my_seq
# Seq('AGTACACTGGT')
# print(my_seq)

import math

from Bio import SeqIO

def matricos_sudarymas(dazniai):
    matrica = [[0.0 for _ in range(8)] for _ in range(8)]
    for n in range(0, 7):
        for m in range(n + 1, 8):
            skaiciavimas = 0.0
            for e in range(len(dazniai[n])):
                skaiciavimas += math.pow((dazniai[n][e] - dazniai[m][e]), 2)
            matrica[n][m] = math.sqrt(skaiciavimas)
            matrica[m][n] = matrica[n][m]
    return matrica


def dikodonu_generavimas():
    aminorugstys = "ARNDCEQGHILKMFPSTWYV"
    dikodonai = []
    for start_kodonas in aminorugstys:
        for pabaigos_kodonas in aminorugstys:
            dikodonai.append(start_kodonas + pabaigos_kodonas)
    # print(dikodonai)
    return dikodonai


def ilgiausia_seka(sekos, stop_kodonas):
    didziausias_ilgis = 0
    didziausia_seka = ""
    for seka in sekos:
        if len(seka) > didziausias_ilgis and str(seka).endswith(stop_kodonas):
            didziausias_ilgis = len(seka)
            didziausia_seka = seka
    return didziausia_seka


def vykdymas(failas):
    for seka in SeqIO.parse("viruses/data/" + failas, "fasta"):
        seka = seka
    # print(failas + " " + seka.seka)

    # 1 uzd
    ilgis = 3
    sekos = [seka.seq[j:j + ilgis] for j in range(len(seka.seq) - 2)]
    sekos_atvirksciai = [seka.seq.reverse_complement()[j:j + ilgis] for j in range(len(seka.seq) - 2)]
    visos_sekos = sekos + sekos_atvirksciai
    # print(visos_sekos)
    k = 0
    kodonu_poros = []
    while k < len(visos_sekos):
        if visos_sekos[k] == 'ATG':
            pradzia = k
            j = k
            while j < len(visos_sekos):
                if visos_sekos[j] == 'TAA' or visos_sekos[j] == 'TAG' \
                        or visos_sekos[j] == 'TGA':
                    pabaiga = j
                    kodonu_poros.append(''.join(str(e) for e in visos_sekos[pradzia:pabaiga + 1]))
                    k = j
                    break
                j += 1
        k += 1
    # print("seka: " + str(seka))

    kodonu_poros = kodonu_poros
    # print(kodonu_poros)

    # 2 uzd
    visos_ilgiausios_sekos = [ilgiausia_seka(kodonu_poros, "TAG"), ilgiausia_seka(kodonu_poros, "TAA"),
                       ilgiausia_seka(kodonu_poros, "TGA")]
    # print(visos_ilgiausios_sekos)

    # 3 uzd
    atfiltruotos_sekos = [seka for seka in kodonu_poros if len(seka) >= 100]
    # print(atfiltruotos_sekos)

    # 4 uzd
    visi_kodonai = ["ATT", "ATC", "ATA", "CTT", "CTC", "CTA", "CTG", "TTA", "TTG", "GTT", "GTC", "GTA", "GTG", "TTT", "TTC",
                  "ATG", "TGT", "TGC", "GCT", "GCC", "GCA", "GCG", "GGT", "GGC", "GGA", "GGG", "CCT", "CCC", "CCA", "CCG",
                  "ACT", "ACC", "ACA", "ACG", "TCT", "TCC", "TCA", "TCG", "AGT", "AGC", "TAT", "TAC", "TGG", "CAA", "CAG",
                  "AAT", "AAC", "CAT", "CAC", "GAA", "GAG", "GAT", "GAC", "AAA", "AAG", "CGT", "CGC", "CGA", "CGG", "AGA",
                  "AGG ,TAA", "TAG", "TGA"]
    visos_sekos = "".join(atfiltruotos_sekos)
    padalinta_seka = [visos_sekos[x:x + 3] for x in range(0, len(visos_sekos), 3)]
    kodonu_dazniai = []
    for kodonas in visi_kodonai:
        skaiciavimas = 0
        for seka in padalinta_seka:
            if seka == kodonas:
                skaiciavimas += 1
        kodonu_dazniai.append(skaiciavimas / len(padalinta_seka))
    kodonu_dazniai = kodonu_dazniai
    # print(kodonu_dazniai)

    from Bio.Seq import Seq
    string = ""
    for el in atfiltruotos_sekos:
        string += el
    dnr = Seq(string)
    isverstas_dnr = dnr.translate()
    # print(isverstas_dnr)
    skaiciav = 0
    dikodonu_seka = []
    while skaiciav < len(isverstas_dnr) - 1:
        if isverstas_dnr[skaiciav] != '*' and isverstas_dnr[skaiciav + 1] != '*':
            dikodonu_seka.append(isverstas_dnr[skaiciav] + isverstas_dnr[skaiciav + 1])
        skaiciav += 1
    # print(dikodonu_seka)
    visi_dikodonai = dikodonu_generavimas()
    dikodonu_dazniai = []
    for dikodonas in visi_dikodonai:
        skaiciavimas = 0
        for seka in dikodonu_seka:
            if seka == dikodonas:
                skaiciavimas += 1
        dikodonu_dazniai.append(skaiciavimas / len(dikodonu_seka))
    dikodonu_dazniai = dikodonu_dazniai
    # print(dikodonu_dazniai)
    return kodonu_dazniai, dikodonu_dazniai


if __name__ == '__main__':
    kodonai, dikodonai = ["ATT", "ATC", "ATA", "CTT", "CTC", "CTA", "CTG", "TTA", "TTG", "GTT", "GTC", "GTA", "GTG", "TTT", "TTC",
                        "ATG", "TGT", "TGC", "GCT", "GCC", "GCA", "GCG", "GGT", "GGC", "GGA", "GGG", "CCT", "CCC", "CCA", "CCG",
                        "ACT", "ACC", "ACA", "ACG", "TCT", "TCC", "TCA", "TCG", "AGT", "AGC", "TAT", "TAC", "TGG", "CAA", "CAG",
                        "AAT", "AAC", "CAT", "CAC", "GAA", "GAG", "GAT", "GAC", "AAA", "AAG", "CGT", "CGC", "CGA", "CGG", "AGA",
                        "AGG ,TAA", "TAG", "TGA"], dikodonu_generavimas()
    visi_dikodonu_dazniai = []
    visi_kodonu_dazniai = []

    for i in range(1, 5):
        kodonai, dikodonai = vykdymas("bacterial" + str(i) + ".fasta")
        visi_kodonu_dazniai.append(kodonai)
        visi_dikodonu_dazniai.append(dikodonai)
    # print(visi_kodonu_dazniai[0][1])
    for i in range(1, 5):
        codonai, dicodonai = vykdymas("mamalian" + str(i) + ".fasta")
        visi_kodonu_dazniai.append(codonai)
        visi_dikodonu_dazniai.append(dicodonai)
    kodonu_matrica = matricos_sudarymas(visi_kodonu_dazniai)
    print("b1 b2 b3 b4 m1 m2 m3 m4")
    for el in kodonu_matrica:
        print(el)
    dikodonu_matrica = matricos_sudarymas(visi_dikodonu_dazniai)
    print("b1 b2 b3 b4 m1 m2 m3 m4")
    for el in dikodonu_matrica:
        print(el)
