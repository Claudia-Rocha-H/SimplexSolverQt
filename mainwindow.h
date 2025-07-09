#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include "modelo.h"
#include "simplexsolver.h"
#include "uibuilder.h"


QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_btnGenerar_clicked();
    void on_btnCalcular_clicked();
    void on_btnSiguiente_clicked();
    void on_btnAnterior_clicked();


private:
    Ui::MainWindow *ui;
    SimplexSolver *simplexSolver;
    UiBuilder *uiBuilder;
    QVector<double> tablaObjetivo;
    QVector<QVector<double>> tablaRestricciones;
    QVector<QString> tablaSignos;
    QVector<double> tablaLD;
};

#endif
