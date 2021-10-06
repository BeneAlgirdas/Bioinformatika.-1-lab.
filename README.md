# Algirdas Benetis
# Bioinformatika. Laboratorinis #1

## Ataskaita

### Užduotis, išpildytos realizacijos, punktai:
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

### Kodonų matrica:
![Kodonų matrica](https://raw.githubusercontent.com/BeneAlgirdas/Bioinformatika.-1-lab./master/images/Kodon%C5%B3%20matrica.png)

### Kodonų matrica - Phylip formatas
![Kodonų matrica - Phylip formatas](https://raw.githubusercontent.com/BeneAlgirdas/Bioinformatika.-1-lab./master/images/Kodon%C5%B3%20matrica%20-%20Phylip%20formatas.png)

### Kodonų atstumų matricos medis
![Kodonų matrica - Phylip formatas](https://raw.githubusercontent.com/BeneAlgirdas/Bioinformatika.-1-lab./master/images/Kodon%C5%B3%20atstum%C5%B3%20matricos%20medis.png)

### Dikodonų matrica:
![Dikodonų matrica](https://raw.githubusercontent.com/BeneAlgirdas/Bioinformatika.-1-lab./master/images/Dikodon%C5%B3%20matrica.png)

### Dikodonų matrica - Phylip formatas
![Dikodonų matrica - Phylip formatas](https://raw.githubusercontent.com/BeneAlgirdas/Bioinformatika.-1-lab./master/images/Dikodon%C5%B3%20matrcia%20-%20Phylip%20formatas.png)

### Dikodonų atstumų matricos medis
![Dikodonų matrica - Phylip formatas](https://raw.githubusercontent.com/BeneAlgirdas/Bioinformatika.-1-lab./master/images/Dikodon%C5%B3%20atstum%C5%B3%20matricos%20medis.png)

#### Palyginimui naudojama formulė: ![Palyginimui naudojama formulė](https://raw.githubusercontent.com/BeneAlgirdas/Bioinformatika.-1-lab./master/images/Palyginimui%20naudojama%20formul%C4%97.png)

