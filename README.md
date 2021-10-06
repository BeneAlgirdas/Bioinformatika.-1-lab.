# Algirdas Benetis
# Bioinformatika.-1-lab.

# Ataskaita

![Kodonų matrica:](https://i7.imageban.ru/out/2021/10/06/dc832591aac782dc9887d285ce5ba17a.png)

# Laboratorinis #1

Instrukcijos
Sveiki,
Užduoties tikslas: įvertinti kodonų ir dikodonų dažnio skirtumus zinduolių ir bakterijų virusuose. Prisegtuke pateikti po keturis  zinduoliu ("mamalian") ir bakteriju virusai (fasta formatas). Parasykite programa, kuri:  
Pateiktoje sekoje fasta formatu surastu visas start ir stop kodonų poras, tarp kurių nebutu stop kodono (ir tiesioginei sekai ir jos reverse komplementui). 
Kiekvienam stop kodonui parinkti toliausiai nuo jo esanti start kodoną (su salyga, kad tarp ju nera kito stop kodono)
Atfiltruokite visus fragmentus ("tai butu baltymų koduojancios sekos"), kurie trumpesni nei 100 fragmentų.
Parasykite funkcijas, kurios ivertintu kodonu ir dikodonu daznius (visi imanomi kodonai/dikodonai ir jų atitinkamas daznis  - gali buti nemazai nuliu, jei ju sekoje nerasite).
Palyginkite kodonu bei dikodonu daznius tarp visu seku (atstumu matrica - kokia formule naudosite/kaip apskaiciuosite - parasykite ataskaitoje).
Ivertinkite, ar bakteriniai ir zinduoliu virusai sudaro atskirus klasterius vertinant kodonu/dikodonu dažniu aspektu. Siuilau atstumu matrica issaugoti tokiu formatu:

5 
Alpha 0.000 1.000 2.000 3.000 3.000 
Beta 1.000 0.000 2.000 3.000 3.000 
Gamma 2.000 2.000 0.000 3.000 3.000 
Delta 3.000 3.000 3.000 0.000 1.000 
Epsilon 3.000 3.000 3.000 1.000 0.000

(tai yra   Phylip formatas)

Pirmas skaičius - klasterizuojamų objektų skaičius), matricos vertės - atstumas tarp objektų.
Tamstų reikalas parinkti atstumo funkciją kodonų bei dikodonų atvejams.

Gavus matricą,  gauti atitinkamą medį, rodantį atitinkamą klasterizavimą neighbour joining metodu,
galite šiame puslapyje:

http://www.trex.uqam.ca/index.php?action=trex&menuD=1&method=2


Taigi iš jusu laukiu:
gitHub repozitorijos su kodu (pasidalinkit - mano email: galzbutas@gmail.com)
Laisvos formos ataskaitos, kurioje:
aprašykite, kaip skaičiavote atstumo funkciją;
kokie medžiai gavosi su kodonais ir dikodonais;
ar skiriasi kodonų ir dikodonų dažnis tarp žinduolių ir bakterijų virusų, kaip klasterizuojasi virusai. Gal kažkuris virusas labai išsiskyrė? Kokie kodonai/dikodonai labiausiai varijuoja?
Sėkmės.
Gediminas

P.S. Nežinau teisingo ataskymo...Įdomu kaip gausis:). Be to, viena iš pateiktų sekų - mielasis covid-19...
