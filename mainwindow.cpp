#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "simplexsolver.h"
#include <QMessageBox>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);
    simplexSolver = new SimplexSolver(ui, this);
    uiBuilder = new UiBuilder(ui, this);
    connect(ui->btnGenerar_2, &QPushButton::clicked, this, &MainWindow::on_btnGenerar_clicked);
    connect(ui->btnCalcular, &QPushButton::clicked, this, &MainWindow::on_btnCalcular_clicked);
    connect(ui->btnSiguiente, &QPushButton::clicked, this, &MainWindow::on_btnSiguiente_clicked);
    connect(ui->btnAnterior_2, &QPushButton::clicked, this, &MainWindow::on_btnAnterior_clicked);
    connect(ui->btnSiguiente, &QPushButton::clicked, this, &MainWindow::on_btnSiguiente_clicked);
}

MainWindow::~MainWindow() {
    delete ui;
}

void MainWindow::on_btnGenerar_clicked() {
    simplexSolver-> reiniciarEstado();
    uiBuilder->generarCampos();

}

void MainWindow::on_btnAnterior_clicked() {
    simplexSolver->pasoAnterior();
}

void MainWindow::on_btnCalcular_clicked() {
    auto datos = uiBuilder->leerDatos();
    if (!datos.ok) return;

    this->tablaObjetivo = datos.objetivo;
    this->tablaRestricciones = datos.restricciones;
    this->tablaSignos = datos.signos;
    this->tablaLD = datos.ladosDerechos;
    simplexSolver->setDatos(datos.objetivo, datos.restricciones, datos.signos, datos.ladosDerechos);
    simplexSolver->mostrarTablaSimplex();
}

void MainWindow::on_btnSiguiente_clicked() {
    simplexSolver->pasoSiguiente();

}


