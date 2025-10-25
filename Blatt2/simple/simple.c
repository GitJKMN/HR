/*
** simple error demonstration to demonstrate power of valgrind
** Julian M. Kunkel - 17.04.2008
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int *mistakes1(void) {
  int buf[] = {1, 1, 2, 3, 4, 5};
  int *buf_ptr = malloc(sizeof(buf)); //Explizite Allokierung von Speicher
  memcpy(buf_ptr, buf, sizeof(buf));  //int-Array buf wird aus dem Stack-Frame in den allokierten Speicher kopiert
  return buf_ptr;
}

int *mistakes2(void) {
  int *buf = malloc(sizeof(int) * 4); //allokierter Speicher muss an die Größe von int angepasst werden
  buf[1] = 2;                         //*p[1] == &mistakes2()[1] ruft zweite Stelle des Arrays auf => buf[1]
  return buf;
}

int *mistakes3(void) {
  /* In dieser Funktion darf kein Speicher direkt d.h. explizit allokiert werden. */
  int mistakes2_ = 0;
  int *buf = mistakes2(); //implizite Allokierung von Speicher über mistakes2
  buf[mistakes2_] = 3;
  return buf;
}

int *mistakes4(void) {
  int *buf = malloc(sizeof(int) * 4); //parallel zu mistakes1()
  buf[0] = 4;                         //*p[3] == mistakes4() ruft erste Stelle des Arrays auf => buf[0]
  //free(buf);
  return buf;
}

int *mistakes5(void) {
  int *buf = malloc(5 * sizeof(int)); //parallel zu mistakes1()
  buf[4] = 5;                         //muss an mistakes5()+4 angepasst werden => buf[4]
  return buf;
}

int main(void) {
  /* Diese Zeile darf NICHT verändert werden! */
  int *p[5] = {&mistakes1()[1], &mistakes2()[1], mistakes3(), mistakes4(), mistakes5()+4};

  /* Die printf-Aufrufe dürfen NICHT verändert werden*/
  printf("1: %d\n", *p[0]);
  printf("2: %d\n", *p[1]);
  printf("3: %d\n", *p[2]);
  printf("4: %d\n", *p[3]);
  printf("5: %d\n", *p[4]);

  /* TODO */
  /* Fügen sie hier die korrekten aufrufe von free() ein */
  //free() muss auf die allgemeine Array-Referenz angewendet werden 
  free(p[0] - 1); 
  free(p[1] - 1);
  free(p[2]);
  free(p[3]);
  free(p[4] - 4);

  return 0;
}
