#ifndef SIMPLEXSOLVER_H
#define SIMPLEXSOLVER_H
#include <QObject>
#include <QVector>
#include <QString>
#include "ui_mainwindow.h"


class SimplexSolver : public QObject {
    Q_OBJECT

public:
    explicit SimplexSolver(Ui::MainWindow *ui, QObject *parent = nullptr);

    void setDatos(const QVector<double>& objetivo,
                  const QVector<QVector<double>>& restricciones,
                  const QVector<QString>& signos,
                  const QVector<double>& ladosDerechos);

    void mostrarTablaSimplex();
    void siguientePaso();
    void actualizarFilaZ();
    void cargarTabla(int paso);
    void pasoSiguiente();
    void pasoAnterior();
    void reiniciarEstado();
    void actualizarLabelPaso();
    void dibujarGrafico();
    void dibujarGrafico3D();
    QVector<QVector3D> calcularVerticesFactibles3D();

private:
    Ui::MainWindow *ui;
    QVector<double> objetivo;
    QVector<QVector<double>> restricciones;
    QVector<QString> signos;
    QVector<double> ladosDerechos;
    int columnaPivote;
    int pasoActual = 0;
    QVector<QVector<QString>> tablasSimplex;
    QWidget* parentWidget;
    void guardarTablaActual();
    void imprimirTablaConsola();
    void mostrarAnalisisSensibilidad();
    bool solucionAlcanzada = false;
    bool esMaximizacion = false;
};

#endif

