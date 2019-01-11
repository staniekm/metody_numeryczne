#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cmath>

using namespace std;

const double eps = 1e-12;//stała przybliżenia zera

bool gauss(int n, double ** AB, double * X)
{
    int i,j,k;
    double m,s;
    
    for(i = 0; i < n - 1; i++)//eliminacja współczynników
    {
        for(j = i + 1; j < n; j++)
        {
            if(fabs(AB[i][i]) < eps) return false;
            m = -AB[j][i] / AB[i][i];
            for(k = i + 1; k <= n; k++)
                AB[j][k] += m * AB[i][k];
        }
    }
    for(i = n - 1; i >= 0; i--)//wyliczanie niewiadomych
    {
        s = AB[i][n];
        for(j = n - 1; j >= i + 1; j--)
            s -= AB[i][j] * X[j];
        if(fabs(AB[i][i]) < eps) return false;
        X[i] = s / AB[i][i];
    }
    return true;
}
int main()
{
    int i, j, k;//zmienne używane przy pętlach
    int z, w = 0;//wybor menu i numer wstrzasu
    int n=4;//ilość niewiadomych
    cout << "podaj który wstrzas chesz obliczyc: ";
    cin >> z;
    switch(z)
    {
        case 1:
            w=0;
            break;
        case 2:
            w=8;
            break;
        case 3:
            w=16;
            break;
        case 4:
            w=24;
            break;
    }
        //zczytanie danych z pliku do macierzy 56x1-------------------------------------------------------------------------------------------------------
    string macierz_danych[57];//stringi z pliku przed konwersja na double
    double macierz_danych_double[57];//double po kowersji
        ifstream dane; //otwarcie pliku z danymi
        dane.open("/Users/stanik/Desktop/studies/sem3/metody_numeryczne/mn_projekt_1_lokalizacja_wstrząsu/mn_projekt_1_lokalizacja_wstrząsu/dane");//ścieżka
        if(dane.is_open())//odczyt danych z pliku do macierzy
        {
            cout << "wczytane dane:" << endl;
            for(i = 0; i < 57; i++)
            {
                getline(dane, macierz_danych[i]);
                cout << macierz_danych[i] << endl; //wyswietlenie z pliku stringow
            }
        }
        for(i = 0; i < 57; i++)//konwersja na stringow na inty
        {
            stringstream s;
            s << macierz_danych[i];
            s >> macierz_danych_double[i];
            //cout << macierz_danych_double[j] << endl;//wyswietlenie czy te same liczby
        }
    
    double suma;//sprawdzenie czy da sie sumowac inty po konwersji czyli czy konwersja działa ;)
        suma=macierz_danych_double[0]+macierz_danych_double[1];
        cout <<"sprawdzenie czy dodaje ('jesli jest 4087 to oki') "<< "'"<< suma << "'" << endl;
    cout << endl;
    
    
    
        //przygotowanie macierzy A-------------------------------------------------------------------------------------------------------------------------
    cout << "Macierz A 7x4" << endl;
    double macierza[7][4];
    for(i=0; i<7; i++)//na wiersze
    {
        double kolumna1=0, kolumna2=0, kolumna3=0, kolumna4=0;
        for(int j=0; j<4; j++)//na kolumny
        {
            if(j==0)
            {
                kolumna1=(-2*(macierz_danych_double[i+1])+2*(macierz_danych_double[0]));//kolumna z x0
                macierza[i][j]=kolumna1;
            }
            if(j==1)
            {
                kolumna2=(-2*(macierz_danych_double[i+9])+2*(macierz_danych_double[8]));// kolumna z y0
                macierza[i][j]=kolumna2;
            }
            if(j==2)
            {
                kolumna3=(-2*(macierz_danych_double[i+17])+2*(macierz_danych_double[16]));// kolumna z z0
                macierza[i][j]=kolumna3;
            }
            if(j==3)
            {
                kolumna4=(pow(macierz_danych_double[56], 2)*(2*macierz_danych_double[i+25+w]-2*macierz_danych_double[24+w])); // kolumna z t0
                macierza[i][j]=kolumna4;
            }
            cout << macierza[i][j] << " ";
        }
        cout << endl;
        //cout << kolumna1 << " " << kolumna2 << " " << kolumna3 << " " << kolumna4 << endl;//obliczone dane(jesli sie powtarzają to poprawnie zapisuje
    }
    cout << endl;
    
    //transponowanie macierzy A----------------------------------------------------------------------------------------------------------------------------
    cout << "Macierz AT 4x7" << endl;
    double macierzat[4][7];
    for(i = 0; i < 7; i++) //transponowanie macierzy
    {
        for(j = 0;j < 4; j++)
        {
            macierzat[j][i]=macierza[i][j];
        }
    }
    for(i = 0; i < 4; i++) //wyswietlanie macierzy transponowanej
    {
        for(j = 0;j < 7; j++)
        {
            cout << macierzat[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    
        //przygotowanie macierzy B-------------------------------------------------------------------------------------------------------------------------
    cout << "Macierz B 7x1" << endl;
    double macierzb[7];//macierz B 7x1
    for(i=0; i<7; i++)//obliczenie i przypisanie pod kolumna1
    {
        double kolumna1;
        kolumna1=-pow(macierz_danych_double[i+1],2)+pow(macierz_danych_double[0],2)-pow(macierz_danych_double[i+9],2)+pow(macierz_danych_double[8],2)-pow(macierz_danych_double[i+17],2)+pow(macierz_danych_double[16],2)+(pow(macierz_danych_double[56],2)*(pow(macierz_danych_double[i+25+w],2)-pow(macierz_danych_double[24+w],2)));
        macierzb[i]=kolumna1;//stworzenie macierzy B z danymi
        cout << macierzb[i]<< endl; //wyswietlenie macierzy b
    }
    cout << endl;
    
        //AT*A--------------------------------------------------------------------------------------------------------------------------
    cout << "AT*A 4x7*7x4=4x4" << endl;
    double macierzata[4][4];
    for(i=0; i<4; i++)//mnożenie
    {
        for(j=0; j<4; j++)
        {
            macierzata[i][j]=0;
            for(k=0; k<7; k++)
            {
                macierzata[i][j] += macierzat[i][k] * macierza[k][j];
            }
        }
    }
    for(i=0; i<4; i++)//wyswietlenie A*AT
    {
        for(j=0; j<4; j++)
        {
            cout << macierzata[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    
        //AT*B-----------------------------------------------------------------------------------------------------------------------
    cout << "AT*B 4x7*7x1=4x1" << endl;
    double macierzatb[7];
    for(i=0; i<4; i++)//mnożenie
    {
        for(j=0; j<1; j++)
        {
            macierzatb[i]=0;
            for(k=0; k<7; k++)
            {
                macierzatb[i]+= macierzat[i][k] * macierzb[k];
            }
        }
    }
    for(i=0; i<4; i++)//wyswietlenie AT*B
    {
            cout << macierzatb[i] << endl;
    }
    cout << endl;
    
        //rozwiązanie metoda elimacji gaussa 4 niewiadome---------------------------------------------------------------------------------
    double **AB, *X;
    
    AB = new double * [n];//tworzymy macierze AB i X
    X  = new double [n];
    
    for(i = 0; i < n; i++)//tworzymy macierz AB dynamiczną
    {
        AB[i] = new double[n + 1];
    }

    for(i = 0; i < 4; i++)//odczytujemy dane dla macierzy AB
    {
        for(j = 0; j < 5; j++)
        {
            if(j<4)
            {
                AB[i][j]=macierzata[i][j];
            }
            else
            {
                AB[i][j]=macierzatb[i];
            }
            
        }
    }
    
    cout << "AB 4x5" << endl;
    for(i = 0; i < 4; i++)//wyswietlenie AB
    {
        for(j = 0; j < 5; j++)
        {
            cout << AB[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    cout << "wyniki: " << endl;
    if(gauss(n,AB,X))//wyliczenie gaussem
    {
        for(i = 0; i < n; i++)
        {
            if(i==0)
            {
                 cout << "x= " << X[i] << endl;//wspolrzedne x
            }
            if(i==1)
            {
                 cout << "y= " << X[i] << endl;//wspolrzedne y
            }
            if(i==2)
            {
                 cout << "z= " << X[i] << endl;//wspolrzedne z
            }
            if(i==3)
            {
                 cout << "t= " << X[i] << endl;//czas wejscia t
            }
        }
    }
    else
    {
        cout << "DZIELNIK ZERO\n";
    }
    
    ofstream wyniki;//zapis wynikow do pliku
    wyniki.open("/Users/stanik/Desktop/studies/sem3/metody_numeryczne/mn_projekt_1_lokalizacja_wstrząsu/mn_projekt_1_lokalizacja_wstrząsu/wyniki.txt");
    if( wyniki.is_open())
    {
        wyniki << "Dane: " << endl;
        for(i = 0; i < 57; i++)//wyswietlenie danych
        {
            wyniki << macierz_danych[i] << endl;
        }
        wyniki << endl;
        wyniki << "Wyniki: " << endl;
        for(i = 0; i < n; i++)//zapisanie wynikow do pliku
        {
            if(i==0)
            {
                wyniki << "x= " << X[i] << endl;//wspolrzedne x
            }
            if(i==1)
            {
                wyniki << "y= " << X[i] << endl;//wspolrzedne y
            }
            if(i==2)
            {
                wyniki << "z= " << X[i] << endl;//wspolrzedne z
            }
            if(i==3)
            {
                wyniki << "t= " << X[i] << endl;//czas wejscia t
            }
        }
    }
    else
    {
        cout <<"Nie udało sie otworzyć pliku :(" << endl;//INFO: otwarcie pliku się nie powiodło
    } //else

    for(i = 0; i < n; i++)//usuwamy macierze z pamięci
    {
        delete [] AB[i];
    }
    delete [] AB;
    delete [] X;

    return 0;
    
    
}

