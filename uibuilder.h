#ifndef UIBUILDER_H
#define UIBUILDER_H

#include <QObject>
#include <QVector>
#include <QLineEdit>
#include <QComboBox>
#include "ui_mainwindow.h"

struct DatosEntrada {
    QVector<double> objetivo;
    QVector<QVector<double>> restricciones;
    QVector<QString> signos;
    QVector<double> ladosDerechos;
    bool ok = false;
};

class UiBuilder : public QObject {
    Q_OBJECT

public:
    explicit UiBuilder(Ui::MainWindow *ui, QObject *parent = nullptr);

    void generarCampos();
    DatosEntrada leerDatos();

private:
    Ui::MainWindow *ui;
    QVector<QLineEdit*> coefObjetivo;
    QVector<QVector<QLineEdit*>> coefRestricciones;
    QVector<QComboBox*> signoRestricciones;
    QVector<QLineEdit*> ladoDerecho;
};

#endif

