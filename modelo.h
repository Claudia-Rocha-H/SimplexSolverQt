#ifndef MODELO_H
#define MODELO_H
#include <QString>
#include <QVector>

struct Restriccion {
    QVector<double> coeficientes;
    QString tipo;
    double terminoIndependiente;
};

class Modelo {
public:
    QString objetivo;
    QVector<double> coeficientesObjetivo;
    QVector<Restriccion> restricciones;

    void limpiar();
};

#endif
