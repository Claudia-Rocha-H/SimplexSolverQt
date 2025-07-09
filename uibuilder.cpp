#include "uibuilder.h"
#include <QMessageBox>
#include <QLabel>    // Asegúrate de incluir QLabel
#include <QLineEdit> // Asegúrate de incluir QLineEdit
#include <QComboBox> // Asegúrate de incluir QComboBox
#include <QLayout>   // Asegúrate de incluir QLayout (para QLayoutItem)

UiBuilder::UiBuilder(Ui::MainWindow *ui, QObject *parent)
    : QObject(parent), ui(ui) {}

void UiBuilder::generarCampos() {
    int numVariables = ui->spinVariables_2->value();
    int numRestricciones = ui->spinRestricciones_2->value();

    QLayoutItem* item;
    while ((item = ui->layoutObjetivo_2->takeAt(0)) != nullptr) {
        if (item->widget()) {
            delete item->widget();
        }
        delete item;
    }
    while ((item = ui->layoutRestricciones_2->takeAt(0)) != nullptr) {
        if (item->widget()) {
            delete item->widget();
        }
        delete item;
    }

    coefObjetivo.clear();
    coefRestricciones.clear();
    signoRestricciones.clear();
    ladoDerecho.clear();

    for (int i = 0; i < numVariables; ++i) {
        QLineEdit* edit = new QLineEdit();
        edit->setFixedWidth(70);
        edit->setPlaceholderText("0");
        coefObjetivo.append(edit);

        ui->layoutObjetivo_2->addWidget(edit, 0, i * 3);

        QLabel* label = new QLabel("x" + QString::number(i + 1));
        label->setContentsMargins(5, 0, 10, 0);
        ui->layoutObjetivo_2->addWidget(label, 0, i * 3 + 1);

        if (i < numVariables - 1) {
            QLabel* plus = new QLabel("+");
            plus->setContentsMargins(0, 0, 10, 0);
            ui->layoutObjetivo_2->addWidget(plus, 0, i * 3 + 2);
        }
    }

    for (int i = 0; i < numRestricciones; ++i) {
        QVector<QLineEdit*> fila;
        for (int j = 0; j < numVariables; ++j) {
            QLineEdit* edit = new QLineEdit();
            edit->setFixedWidth(70);
            edit->setPlaceholderText("0");
            fila.append(edit);
            ui->layoutRestricciones_2->addWidget(edit, i, j * 3);

            QLabel* label = new QLabel("x" + QString::number(j + 1));
            label->setContentsMargins(5, 0, 10, 0);
            ui->layoutRestricciones_2->addWidget(label, i, j * 3 + 1);

            if (j < numVariables - 1) {
                QLabel* plus = new QLabel("+");
                plus->setContentsMargins(0, 0, 10, 0);
                ui->layoutRestricciones_2->addWidget(plus, i, j * 3 + 2);
            }
        }

        coefRestricciones.append(fila);

        QComboBox* signo = new QComboBox();
        signo->addItems({"<=", "=", ">="});
        signo->setFixedWidth(80);
        signoRestricciones.append(signo);
        ui->layoutRestricciones_2->addWidget(signo, i, numVariables * 3);

        QLineEdit* ld = new QLineEdit();
        ld->setFixedWidth(80);
        ld->setPlaceholderText("0");
        ladoDerecho.append(ld);
        ui->layoutRestricciones_2->addWidget(ld, i, numVariables * 3 + 1);
    }
}

DatosEntrada UiBuilder::leerDatos() {
    DatosEntrada datos;
    bool todoOk = true;

    for (auto* edit : coefObjetivo) {
        QString texto = edit->text().trimmed();
        if (texto.isEmpty()) {
            QMessageBox::warning(nullptr, "Error", "Por favor llena todos los campos de la función objetivo.");
            todoOk = false;
            break;
        }
        bool ok;
        double valor = texto.toDouble(&ok);
        if (!ok) {
            QMessageBox::warning(nullptr, "Error", "Coeficiente inválido en la función objetivo: " + texto);
            todoOk = false;
            break;
        }
        datos.objetivo.append(valor);
    }

    for (int i = 0; i < coefRestricciones.size() && todoOk; ++i) {
        QVector<double> fila;
        for (auto* edit : coefRestricciones[i]) {
            QString texto = edit->text().trimmed();
            if (texto.isEmpty()) {
                fila.append(0.0);
            } else {
                bool ok;
                double valor = texto.toDouble(&ok);
                if (!ok) {
                    QMessageBox::warning(nullptr, "Error", "Coeficiente inválido en restricción " + QString::number(i + 1) + ": " + texto);
                    todoOk = false;
                    break;
                }
                fila.append(valor);
            }
        }
        datos.restricciones.append(fila);
        datos.signos.append(signoRestricciones[i]->currentText());

        QString textoLD = ladoDerecho[i]->text().trimmed();
        if (textoLD.isEmpty()) {
            QMessageBox::warning(nullptr, "Error", "Lado derecho vacío en restricción " + QString::number(i + 1));
            todoOk = false;
            break;
        }
        bool ok;
        double valorLD = textoLD.toDouble(&ok);
        if (!ok) {
            QMessageBox::warning(nullptr, "Error", "Lado derecho inválido en restricción " + QString::number(i + 1) + ": " + textoLD);
            todoOk = false;
            break;
        }
        datos.ladosDerechos.append(valorLD);
    }

    datos.ok = todoOk;
    return datos;
}
