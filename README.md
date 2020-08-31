Mój program przedstawia rozwiązanie niestacjonarne dla schematu
dwupunktowego 2D. Znajduję się w nim sześć klas: Element, Node,
Grid, GlobalData, Matrix i Main. Moje oprogramowanie zaczęłam od
wygenerowania siatki (dla pierwszego testCase siatka 9-elementowa) z
numeracją od zera dla elementów i węzłów.

Każdy element siatki zawiera dwie 4-elementowe tablice: ID (id każdego
z 4 węzłów elementu), BC (flagę każdego węzła, która przyjmuje
wartości true lub false) oraz swoje własne id (ownID).
Każdy węzeł zawiera swoje współrzędne (x, y), temperature (t), status
węzła (BC) oraz swoje id (ownID).
Klasy Element i Node zawierają konstruktory oraz kilka getter’,ów i
setter’ów. Dodatkowo w klasie Elemnnt tworze metody wyświetlania
(show, printHbc, printLocalC, printLocalH).
W klasie GlobalData przechowuje wszystkie dane: początkową
temperaturę (initialT), alfę (alfa), wysokość oraz szerokość siatki (H, W),
ilość węzłów na wysokość oraz szerokość (nH, nW), ciepło właściwe
(heat), gęstość (density), współczynnik przewodzenia ciepla (k),
temperature otoczenia (ambientT), czas kroku symulacji (stepTime),
czas symulacji (simulationTime) oraz ilość wszystkich węzłów i
elementów (numberOfN, numberOfE), które otrzymuje ze wzrów:
● nH * nW
● ( nH - 1) * ( nW - 1)
W klasie Matrix tworzę konstruktor oraz metody:
● MatrixPcx() - która zawiera tablicę Pcx[4][2] z wartościami moich
czterech punktów całkowania
● MatrixN_Ksi() - która przyjmuje wartość eta i zawiera 4
elementową tablicę N_Ksi(), ze wzorami funkcji kształtu po ksi
● MatrixN_Eta() - która przyjmuje wartość ksi i zawiera 4
elementową tablicę N_Eta, ze wzorami funkcji kształtu po eta
● MatrixNx() - która przyjmuje wartość ksi oraz eta i zawiera 4
elementową tablicę ze wzorami funkcji kształtu
● generateMatrices() - w tej metodzie tworzę sobie
➔ dwuwymiarowe tablice 4x4 (dNdKsi, dNdEta, N), które liczą
mi pochodne funkcji kształtu względem ksi i eta dla
wszystkich 4 punktów całkowania
➔ pętle po elementach, w której pobieram wartości x i y dla
poszczególnych węzłów, w tej pętli tworzę kolejną pętle po
punktach całkowania, w której:
➢ liczę Jakobian
➢ wyznacznik macierzy Jakobiego
➢ pochodne funkcji kształtu po x i y
➢ wyliczamy macierz H (localH)
➢ wyliczamy macierz C (localC),
➢ lokalne macierze H i C (dla każdego punktu
całkowania) sumuje dla każdego elementu
➔ po wyjściu z pętli po punktach całkowania (cały czas
jesteśmy w pętli po elementach) sprawdzam dla każdej
ściany czy zachodzi warunek brzegowy (tworzę w tym celu
dla każdej ściany if), jeśli dla 2 węzłów na ścianie BC = true
(isBC()), nakładamy na tą ścianę warunek brzegowy i
liczymy:
➢ wartości wszystkich czterech funkcji kształtu dla
wartości ksi i eta
➢ detJ - stosunek długości układu globalnego do długości
układu lokalnego (długość układu lokalnego -> 1-(-1) =
2)
➢ wektor P (VP)
➢ macierz Hbc (Hbc)
➔ następnie tworzę agregację macierzy H (globalH), C
(globalC), Hbc (globalHbc) i wektora P (GP)
➔ sumuje globalną macierz H z macierzą Hbc (macierz H to
całka po objętości, a Hbc to całka po powierzchni, żeby
odnieść się do całości musimy to połączyć)
● nextStep() - w której tworzę pętle po czasie, w której to
wyznaczam wektor {t1} poprzez rozwiązanie układu równań (z
każdym krokiem czasowym (delta tau = const) temperatura
początkowa {t0}, będzie przyjmowała wartości {t1})
● lslove() - przyjmuję dwuwymiarową macierz i wektor, zawiera
algorytm rozwiązywania układów równań metodą Gaussa. Metodę
tą wykorzystujemy w metodzie nextStep aby wyznaczyć {t1} w
naszym równaniu:
[H] * {t1} = {P}
zwraca wektor {t1}
● getMinTemperature() - przechodzi po elementach w tablicy i
zwraca minimum
● getMaxTemperature() - przechodzi po elementach w tablicy i
zwraca maximum
