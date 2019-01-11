//
//  main.cpp
//  mn_projekt_2_parametry_drgań
//
//  Created by Marcin Staniek on 19/11/2018.
//  Copyright © 2018 Marcin Staniek. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cmath>

using namespace std;

const double eps = 1e-12;//stała przybliżenia zera

bool gauss(int n, double ** GAUSS, double * B)
{
    int i,j,k;
    double m,s;
    
    for(i = 0; i < n - 1; i++)//eliminacja współczynników
    {
        for(j = i + 1; j < n; j++)
        {
            if(fabs(GAUSS[i][i]) < eps) return false;
            m = -GAUSS[j][i] / GAUSS[i][i];
            for(k = i + 1; k <= n; k++)
                GAUSS[j][k] += m * GAUSS[i][k];
        }
    }
    for(i = n - 1; i >= 0; i--)//wyliczanie niewiadomych
    {
        s = GAUSS[i][n];
        for(j = n - 1; j >= i + 1; j--)
            s -= GAUSS[i][j] * B[j];
        if(fabs(GAUSS[i][i]) < eps) return false;
        B[i] = s / GAUSS[i][i];
    }
    return true;
}
int main()
{
    int n=4;//ilosc niewiadomych do gaussa
    int i, j, l;//zmienne w petlach
    int w=0, k=0;//ilosci wierszy i kolumn
    //double energia_sejsmiczna_wstrzasu=200000;
    double macierz_wspolrzednych_x[3], macierz_wspolrzednych_y[3];
    macierz_wspolrzednych_x[0]=3728;//wspolrzedne x stanowisk 1  2 3
    macierz_wspolrzednych_x[1]=3728;
    macierz_wspolrzednych_x[2]=640;
    
    macierz_wspolrzednych_y[0]=426;//wpolrzedne y stanowisk 1 2 3
    macierz_wspolrzednych_y[1]=2494;
    macierz_wspolrzednych_y[2]=2688;
    
    //zczytanie danych z pliku do macierzy 60x7------------------------------------------------------------------------------------------------
    cout << "wczytane dane:" << endl;
    w=60;//
    k=7;
    string macierz_danych[w][k];//stringi z pliku przed konwersja na double
    double macierz_danych_double[w][k];//double po kowersji
    ifstream dane; //otwarcie pliku z danymi
    dane.open("/Users/stanik/Desktop/studies/sem3/metody_numeryczne/mn_projekt_2_parametry_drgań/mn_projekt_2_parametry_drgań/dane");//ścieżka
    if(dane.is_open())//odczyt danych z pliku do macierzy
    {
        for(i = 0; i < w; i++)
        {
            for(j = 0; j < k; j++)
            {
                dane >> macierz_danych[i][j];
                cout << macierz_danych[i][j] << " ";
            }
            cout << endl;
        }
    }
    cout << endl;
    for(i = 0; i < w; i++)//konwersja na stringow na liczby
    {
        for(j = 0; j < k; j++)
        {
            stringstream s;
            s << macierz_danych[i][j];
            s >> macierz_danych_double[i][j];
            cout << macierz_danych_double[i][j] << " "; //wyswietlenie czy te same liczby
        }
        cout << endl;
    }
    cout << endl;
    double suma;//sprawdzenie czy da sie sumowac inty po konwersji czyli czy konwersja działa ;)
    suma=macierz_danych_double[0][0]+macierz_danych_double[1][1];
    cout <<"sprawdzenie czy dodaje ('jesli jest 3 to oki') "<< "'"<< suma << "'" << endl;
    cout << endl;
    
    //A 60x4----------------------------------------------------------------------------------------------------------------
    cout << "A 60x4: " << endl;
    w=60;
    k=4;
    double macierza[w][k];
    for (i=0; i<w; i++)//tworzenie macierzy A
    {
        for (j=0; j<k; j++)
        {
            macierza[i][j]=0;
            if(j==0)
            {
                macierza[i][j]=log10(macierz_danych_double[i][2]);
            }
            if(j==1)//kolumna druga zmienia sie co 20 bo 3 wstrzasy !!dziwnie zapisane R=log10 powinno byc log10(R)
            {
                if((i >= 0) && (i <= 19))
                {
                macierza[i][j]=log10(sqrt(pow(macierz_danych_double[i][3]-macierz_wspolrzednych_x[0],2) + pow(macierz_danych_double[i][4]-macierz_wspolrzednych_y[0],2) + pow(500,2)));
                }
                if((i >= 20) && (i <= 39))
                {
                macierza[i][j]=log10(sqrt(pow(macierz_danych_double[i][3]-macierz_wspolrzednych_x[1],2) + pow(macierz_danych_double[i][4]-macierz_wspolrzednych_y[1],2) + pow(500,2)));
                }
                if((i >= 40) && (i <= 59))
                {
                macierza[i][j]=log10(sqrt(pow(macierz_danych_double[i][3]-macierz_wspolrzednych_x[2],2) + pow(macierz_danych_double[i][4]-macierz_wspolrzednych_y[2],2) + pow(500,2)));
                }
            }
            if(j==2)//kolumna trzecia zmienia sie co 20 bo 3 wstrzasy
            {
                if((i >= 0) && (i <= 19))
                {
                    macierza[i][j]=sqrt(pow(macierz_danych_double[i][3]-macierz_wspolrzednych_x[0],2) + pow(macierz_danych_double[i][4]-macierz_wspolrzednych_y[0],2) + pow(500,2));
                }
                if((i >= 20) && (i <= 39))
                {
                    macierza[i][j]=sqrt(pow(macierz_danych_double[i][3]-macierz_wspolrzednych_x[1],2) + pow(macierz_danych_double[i][4]-macierz_wspolrzednych_y[1],2) + pow(500,2));
                }
                if((i >= 40) && (i <= 59))
                {
                    macierza[i][j]=sqrt(pow(macierz_danych_double[i][3]-macierz_wspolrzednych_x[2],2) + pow(macierz_danych_double[i][4]-macierz_wspolrzednych_y[2],2) + pow(500,2));
                }
            } 
            if(j==3)
            {
                macierza[i][j]=1;
            }
        }
    }
    for (i=0; i<w; i++)//wyswietlenie macierzy A
    {
        for (j=0; j<k; j++)
        {
            cout << macierza[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    
        //AT 4x60-----------------------------------------------------------------------------------------------------------
    cout << "AT 4x60: " << endl;
    w=4;
    k=60;
    double macierzat[w][k];
    for(i = 0; i < k; i++) //transponowanie macierzy
    {
        for(j = 0;j < w; j++)
        {
            macierzat[j][i]=macierza[i][j];
        }
    }
    for(i = 0; i < w; i++) //wyswietlanie macierzy transponowanej
    {
        for(j = 0;j < k; j++)
        {
            cout << macierzat[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

        //Y 60x1-------------------------------------------------------------------------------------------------------
    cout << "Y 60x1: " << endl;
    w=60;
    k=1;
    double macierzY[w][k];
    for (i=0; i<w; i++)//tworzenie macierzy B
    {
        for (j=0; j<k; j++)
        {
            macierzY[i][j]=0;
            macierzY[i][j]=log10(macierz_danych_double[i][5]);
            cout << macierzY[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    
    //AT*A 4x60*60x4=4x4---------------------------------------------------------------------------------------------------------
    cout << "AT*A 4x60*60x4=4x4:" << endl;
    w=4;
    k=4;
    double macierzata[w][k];
    for(i=0; i<w; i++)//mnożenie
    {
        for(j=0; j<k; j++)
        {
            macierzata[i][j]=0;
            for(l=0; l<60; l++)
            {
                macierzata[i][j] += macierzat[i][l] * macierza[l][j];
            }
        }
    }
    for(i=0; i<w; i++)//wyswietlenie AT*A
    {
        for(j=0; j<k; j++)
        {
            cout << macierzata[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    
    //AT*Y 4X60*60x1=4x1---------------------------------------------------------------------------------------------------------
    cout << "AT*Y 4x60*60x1=4x1" << endl;
    w=4;
    k=1;
    double macierzatY[w];
    for(i=0; i<w; i++)//mnożenie
    {
        for(j=0; j<k; j++)
        {
            macierzatY[i]=0;
            for(l=0; l<60; l++)
            {
                macierzatY[i]+= macierzat[i][l] * macierzY[l][j];
            }
        }
    }
    for(i=0; i<w; i++)//wyswietlenie AT*Y
    {
        for(j=0; j<k ;j++)
        {
            cout << macierzatY[i]<< " ";
        }
        cout << endl;
    }
    cout << endl;
    
    //rozwiazanie gaussem--------------------------------------------------------------------------------------------------------
    double **GAUSS, *B;
    
    GAUSS = new double * [n];//tworzymy macierze GAUSS i B
    B  = new double [n];
    
    for(i = 0; i < n; i++)//tworzymy macierz GAUSS dynamiczną
    {
        GAUSS[i] = new double[n + 1];
    }
    cout << "GAUSS 4x5:" << endl;
    for(i = 0; i < 4; i++)//odczytujemy dane dla macierzy GAUSS
    {
        for(j = 0; j < 5; j++)
        {
            if(j<4)
            {
                GAUSS[i][j]=macierzata[i][j];
            }
            else
            {
                GAUSS[i][j]=macierzatY[i];
            }
            
        }
    }
    
    cout << "Wstawienie do Gauss'a 4x5" << endl;
    for(i = 0; i < 4; i++)//wyswietlenie GAUSS
    {
        for(j = 0; j < 5; j++)
        {
            cout << GAUSS[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    cout << "Rozwiazanie GAUSSA: " << endl;
    if(gauss(n,GAUSS,B))//wyliczenie gaussem
    {
        for(i = 0; i < n; i++)
        {
            if(i==0)
            {
                cout << "b1= " << B[i] << endl;
            }
            if(i==1)
            {
                cout << "b2= " << B[i] << endl;
            }
            if(i==2)
            {
                cout << "b3= " << B[i] << endl;
            }
            if(i==3)
            {
                cout << "b4= " << B[i] << endl;
            }
        }
    }
    else
    {
        cout << "DZIELNIK ZERO\n";
    }
    cout << endl;
    
    double odleglosci_epicentralne[4];
    odleglosci_epicentralne[0]=0;//strefa epicentralna
    odleglosci_epicentralne[1]=500;//odleglosci epicentralne
    odleglosci_epicentralne[2]=1000;
    odleglosci_epicentralne[3]=1500;
    
    double Repicentralne,R500,R1000,R1500;
    Repicentralne=sqrt(pow(odleglosci_epicentralne[0],2)+pow(500, 2));
    R500=sqrt(pow(odleglosci_epicentralne[1],2)+pow(500, 2));
    R1000=sqrt(pow(odleglosci_epicentralne[2],2)+pow(500, 2));
    R1500=sqrt(pow(odleglosci_epicentralne[3],2)+pow(500, 2));
    
    double Aepicentralne,Aepicentralne500, Aepicentralne1000, Aepicentralne1500;
    Aepicentralne=pow(10,(B[0]*log10(200000))+(B[1]*log10(Repicentralne))+(B[2]*Repicentralne)+B[3]);
    Aepicentralne500=pow(10,(B[0]*log10(200000))+(B[1]*log10(R500))+(B[2]*R500)+B[3]);
    Aepicentralne1000=pow(10,(B[0]*log10(200000))+(B[1]*log10(R1000))+(B[2]*R1000)+B[3]);
    Aepicentralne1500=pow(10,(B[0]*log10(200000))+(B[1]*log10(R1500))+(B[2]*R1500)+B[3]);
    
    
    //przedziały ufnosci----------------------------------------------------------------------------------------------------------
    double odchylenie_standardowe[60], odchylenie_standardoweSR=0;//E
    double sigma=0;//o-
    double lambda=1.6707;//h

    //Yprog maciedz 60x1 przyspieszen prognozowanych
    cout << "Yprog 60x1:" << endl;
    w=60;
    k=1;
    double macierzYprog[w][k];//logAuf
    for (i=0; i<60; i++)
    {
        for (j=0; j<1; j++)
        {
            macierzYprog[i][j]=(B[0]*macierza[i][0])+(B[1]*macierza[i][1])+(B[2]*macierza[i][2])+B[3];
            cout << macierzYprog[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    
    //odchylenie standardowe
    cout << "Odchylenie standardowe: " << endl;
    for(i=0;i<60;i++)
    {
        odchylenie_standardowe[i]=macierzY[i][0]+macierzYprog[i][0];
        cout << odchylenie_standardowe[i] << endl;
    }
    cout << endl;
    
    for(i=0;i<60;i++)//odchylenie standardowe średnie
    {
        odchylenie_standardoweSR+=odchylenie_standardowe[i];
    }
    odchylenie_standardoweSR=odchylenie_standardoweSR/60;
    cout << "Odchylenie standardowe średnie: " << endl << odchylenie_standardoweSR << endl;
    cout << endl;
    
    //sigma
    cout << "Sigma: " << endl;
    for(i=0;i<60;i++)
    {
        sigma+=pow(odchylenie_standardowe[i]-odchylenie_standardoweSR,2);
    }
    sigma=sqrt(sigma)/59;
    cout << sigma << endl;
    cout << endl;
    
    double AepicentralneUF,Aepicentralne500UF, Aepicentralne1000UF, Aepicentralne1500UF;
    AepicentralneUF=pow(10,log10(Aepicentralne)+(sigma*lambda));
    Aepicentralne500UF=pow(10,log10(Aepicentralne500)+(sigma*lambda));
    Aepicentralne1000UF=pow(10,log10(Aepicentralne1000)+(sigma*lambda));
    Aepicentralne1500UF=pow(10,log10(Aepicentralne1500)+(sigma*lambda));
    
    cout << "A: " << Aepicentralne << endl;
    cout << "A500: " << Aepicentralne500 <<  endl;
    cout << "A1000: " << Aepicentralne1000 << endl;
    cout << "A1500: " << Aepicentralne1500 << endl;
    cout << endl;
    cout << "Auf: " << AepicentralneUF<< endl;
    cout << "Auf500: " << Aepicentralne500UF <<  endl;
    cout << "Auf1000: " << Aepicentralne1000UF << endl;
    cout << "Auf1500: " << Aepicentralne1500UF << endl;
    
}

