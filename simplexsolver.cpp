#include "simplexsolver.h"
#include <QMessageBox>
#include <QDebug>
#include <limits>
#include <QGraphicsTextItem>
#include <Eigen/Dense>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QFile>
#include <QDir>
#include <QFileInfo>

using namespace Eigen;

SimplexSolver::SimplexSolver(Ui::MainWindow *ui, QObject *parent)
    : QObject(parent), ui(ui) {
    parentWidget = qobject_cast<QWidget*>(parent);
}

void SimplexSolver::setDatos(const QVector<double>& obj,
                             const QVector<QVector<double>>& restr,
                             const QVector<QString>& sig,
                             const QVector<double>& ld) {
    objetivo = obj;
    restricciones = restr;
    signos = sig;
    ladosDerechos = ld;
}

void SimplexSolver::mostrarTablaSimplex() {
    int numVars = objetivo.size();
    int numRestricciones = restricciones.size();

    int numHolguras = 0;
    int numArtificiales = 0;

    for (const QString& signo : signos) {
        if (signo == "<=") numHolguras++;
        else if (signo == ">=" || signo == "=") numArtificiales++;
    }

    int totalColumnas = 3 + numVars + numHolguras + numArtificiales + 2;
    int totalFilas = numRestricciones + 1;

    ui->tablaSimplex->setRowCount(totalFilas);
    ui->tablaSimplex->setColumnCount(totalColumnas);

    QStringList headers;
    headers << "Ec #" << "VB" << "Z";
    for (int i = 0; i < numVars; ++i)
        headers << "x" + QString::number(i + 1);
    for (int i = 0; i < numHolguras; ++i)
        headers << "s" + QString::number(i + 1);
    for (int i = 0; i < numArtificiales; ++i)
        headers << "a" + QString::number(i + 1);
    headers << "LD" << "Razón";

    ui->tablaSimplex->setHorizontalHeaderLabels(headers);

    ui->tablaSimplex->setItem(0, 0, new QTableWidgetItem("(0)"));
    ui->tablaSimplex->setItem(0, 1, new QTableWidgetItem("Z"));
    ui->tablaSimplex->setItem(0, 2, new QTableWidgetItem("1"));
    esMaximizacion = ui->rad_max->isChecked();
    for (int i = 0; i < numVars; ++i) {
        double coef = esMaximizacion ? -objetivo[i] : objetivo[i];
        ui->tablaSimplex->setItem(0, 3 + i, new QTableWidgetItem(QString::number(coef)));
    }

    bool esMaximizacion = ui->rad_max->isChecked();
    double M = 1e6;

    for (int i = 0; i < numHolguras; ++i) {
        ui->tablaSimplex->setItem(0, 3 + numVars + i, new QTableWidgetItem("0"));
    }

    for (int i = 0; i < numArtificiales; ++i) {
        double valorM = esMaximizacion ? -M : M;
        ui->tablaSimplex->setItem(0, 3 + numVars + numHolguras + i, new QTableWidgetItem(QString::number(valorM)));
    }


    ui->tablaSimplex->setItem(0, totalColumnas - 2, new QTableWidgetItem("0"));
    ui->tablaSimplex->setItem(0, totalColumnas - 1, new QTableWidgetItem(""));

    int colPivote = 3;
    double minZ = 0;
    for (int c = 0; c < numVars; ++c) {
        double zCoef = ui->tablaSimplex->item(0, 3 + c)->text().toDouble();
        if (zCoef < minZ) {
            minZ = zCoef;
            colPivote = 3 + c;
        }
    }

    int colHolgura = 0;
    int colArtificial = 0;

    for (int i = 0; i < numRestricciones; ++i) {
        int row = i + 1;

        ui->tablaSimplex->setItem(row, 0, new QTableWidgetItem("(" + QString::number(row) + ")"));

        QString signo = signos[i];
        QString vb;

        if (signo == "<=") {
            vb = "s" + QString::number(colHolgura + 1);
        } else if (signo == ">=" || signo == "=") {
            vb = "a" + QString::number(colArtificial + 1);
        }

        ui->tablaSimplex->setItem(row, 1, new QTableWidgetItem(vb));

        ui->tablaSimplex->setItem(row, 2, new QTableWidgetItem("0"));

        for (int j = 0; j < numVars; ++j) {
            ui->tablaSimplex->setItem(row, 3 + j, new QTableWidgetItem(QString::number(restricciones[i][j])));
        }

        for (int j = 0; j < numHolguras + numArtificiales; ++j) {
            ui->tablaSimplex->setItem(row, 3 + numVars + j, new QTableWidgetItem("0"));
        }

        if (signo == "<=") {
            ui->tablaSimplex->setItem(row, 3 + numVars + colHolgura, new QTableWidgetItem("1"));
            colHolgura++;
        } else if (signo == ">=") {
            ui->tablaSimplex->setItem(row, 3 + numVars + colHolgura, new QTableWidgetItem("-1"));
            ui->tablaSimplex->setItem(row, 3 + numVars + numHolguras + colArtificial, new QTableWidgetItem("1"));
            colHolgura++;
            colArtificial++;
        } else if (signo == "=") {
            ui->tablaSimplex->setItem(row, 3 + numVars + numHolguras + colArtificial, new QTableWidgetItem("1"));
            colArtificial++;
        }

        double LD = ladosDerechos[i];
        ui->tablaSimplex->setItem(row, totalColumnas - 2, new QTableWidgetItem(QString::number(LD)));


        double valorPivote = ui->tablaSimplex->item(row, colPivote)->text().toDouble();
        if (valorPivote > 0)
            ui->tablaSimplex->setItem(row, totalColumnas - 1,
                                      new QTableWidgetItem(QString::number(LD / valorPivote, 'f', 2)));
        else
            ui->tablaSimplex->setItem(row, totalColumnas - 1, new QTableWidgetItem("-"));
    }

    tablasSimplex.clear();
    guardarTablaActual();
    pasoActual = 0;
    actualizarLabelPaso();
    imprimirTablaConsola();
    this->columnaPivote = colPivote;
}

void SimplexSolver::siguientePaso() {
    int filas = ui->tablaSimplex->rowCount();
    int colLD = ui->tablaSimplex->columnCount() - 2;
    int colRazon = ui->tablaSimplex->columnCount() - 1;

    esMaximizacion = ui->rad_max->isChecked();
    bool optimo = true;
    const double epsilon = 1e-6;

    for (int c = 3; c < ui->tablaSimplex->columnCount() - 2; ++c) {
        QString nombreVariable = ui->tablaSimplex->horizontalHeaderItem(c)->text();
        if (nombreVariable.startsWith("x")) {
            double zCoef = ui->tablaSimplex->item(0, c)->text().toDouble();
            if ((esMaximizacion && zCoef < -epsilon) || (!esMaximizacion && zCoef > epsilon)) {
                optimo = false;
                break;
            }
        }
    }

    if (optimo) {
        if (!solucionAlcanzada) {
            guardarTablaActual();
            solucionAlcanzada = true;
        }

        QString mensaje = " Se ha alcanzado la solución óptima.\n\n";
        mensaje += "Valor de Z = " + ui->tablaSimplex->item(0, ui->tablaSimplex->columnCount() - 2)->text() + "\n\n";
        mensaje += "Variables básicas:\n";
        for (int i = 1; i < ui->tablaSimplex->rowCount(); ++i) {
            QString vb = ui->tablaSimplex->item(i, 1)->text();
            QString valor = ui->tablaSimplex->item(i, ui->tablaSimplex->columnCount() - 2)->text();
            mensaje += "   - " + vb + " = " + valor + "\n";
        }

        QMessageBox::information(nullptr, "Resultado óptimo", mensaje);
        if (ui->textEditResultado) {
            ui->textEditResultado->setPlainText(mensaje);
        }
        mostrarAnalisisSensibilidad();

        if (objetivo.size() == 2)
            dibujarGrafico();
        else if (objetivo.size() == 3)
            dibujarGrafico3D();

        return;
    }
    int nuevaColPivote = -1;
    double menorValor = 0.0;
    double mejorValor = esMaximizacion ? std::numeric_limits<double>::max() : -std::numeric_limits<double>::max();

    for (int c = 3; c < ui->tablaSimplex->columnCount() - 2; ++c) {
        QString nombreVariable = ui->tablaSimplex->horizontalHeaderItem(c)->text();
        if (nombreVariable.startsWith("x")) {
            double zCoef = ui->tablaSimplex->item(0, c)->text().toDouble();

            if (esMaximizacion) {
                if (zCoef < mejorValor - epsilon) {
                    mejorValor = zCoef;
                    nuevaColPivote = c;
                }
            } else {
                if (zCoef > mejorValor + epsilon) {
                    mejorValor = zCoef;
                    nuevaColPivote = c;
                }
            }

        }
    }


    if (nuevaColPivote == -1) {
        bool todasArtificialesEliminadas = true;

        for (int i = 1; i < filas; ++i) {
            QString vb = ui->tablaSimplex->item(i, 1)->text();
            if (vb.startsWith("a")) {
                double val = ui->tablaSimplex->item(i, colLD)->text().toDouble();
                if (val > 1e-6) {
                    todasArtificialesEliminadas = false;
                    break;
                }
            }
        }

        if (!todasArtificialesEliminadas && pasoActual > 0) {
            QMessageBox::warning(nullptr, "Problema inviable", "El problema no tiene solución factible: variables artificiales siguen en la base.");
            return;
        }


        if (!solucionAlcanzada) {
            guardarTablaActual();
            solucionAlcanzada = true;
        }

        QString mensaje = " Se ha alcanzado la solución óptima.\n\n";
        mensaje += "Valor de Z = " + ui->tablaSimplex->item(0, colLD)->text() + "\n\n";
        mensaje += "Variables básicas:\n";

        for (int i = 1; i < filas; ++i) {
            QString vb = ui->tablaSimplex->item(i, 1)->text();
            QString valor = ui->tablaSimplex->item(i, colLD)->text();
            mensaje += "   - " + vb + " = " + valor + "\n";
        }

        QMessageBox::information(nullptr, "Resultado óptimo", mensaje);
        if (ui->textEditResultado) {
            ui->textEditResultado->setPlainText(mensaje);
        }
        mostrarAnalisisSensibilidad();

        if (objetivo.size() == 2)
            dibujarGrafico();
        else if (objetivo.size() == 3)
            dibujarGrafico3D();

        return;

    }

    solucionAlcanzada = false;
    this->columnaPivote = nuevaColPivote;

    double menorRazon = std::numeric_limits<double>::max();
    int filaPivote = -1;
    for (int i = 1; i < filas; ++i) {
        double valorColPivote = ui->tablaSimplex->item(i, columnaPivote)->text().toDouble();
        double valorLD = ui->tablaSimplex->item(i, colLD)->text().toDouble();

        if (valorColPivote > 0) {
            double razon = valorLD / valorColPivote;
            if (razon < menorRazon) {
                menorRazon = razon;
                filaPivote = i;
            }
        }
    }

    if (filaPivote == -1) {
        QMessageBox::information(nullptr, "Fin", "No hay más pasos posibles (óptimo alcanzado o sin solución).");
        return;
    }


    double valorPivote = ui->tablaSimplex->item(filaPivote, columnaPivote)->text().toDouble();

    if (qFuzzyCompare(valorPivote, 0.0)) {
        QMessageBox::warning(nullptr, "Error", "Valor pivote es cero, no se puede continuar.");
        return;
    }

    for (int j = 2; j < ui->tablaSimplex->columnCount() - 1; ++j) {
        double valor = ui->tablaSimplex->item(filaPivote, j)->text().toDouble();
        double nuevoValor = valor / valorPivote;
        ui->tablaSimplex->item(filaPivote, j)->setText(QString::number(nuevoValor, 'f', 6));
    }

    for (int i = 0; i < filas; ++i) {
        if (i == filaPivote) continue;

        double factor = ui->tablaSimplex->item(i, columnaPivote)->text().toDouble();

        for (int j = 2; j < ui->tablaSimplex->columnCount() - 1; ++j) {
            double valorActual = ui->tablaSimplex->item(i, j)->text().toDouble();
            double valorPivoteFila = ui->tablaSimplex->item(filaPivote, j)->text().toDouble();
            double nuevoValor = valorActual - factor * valorPivoteFila;
            ui->tablaSimplex->item(i, j)->setText(QString::number(nuevoValor, 'f', 6));
        }
    }

    QString encabezado = ui->tablaSimplex->horizontalHeaderItem(columnaPivote)->text();
    ui->tablaSimplex->setItem(filaPivote, 1, new QTableWidgetItem(encabezado));

    for (int i = 1; i < filas; ++i) {
        double valorColPivote = ui->tablaSimplex->item(i, columnaPivote)->text().toDouble();
        double valorLD = ui->tablaSimplex->item(i, colLD)->text().toDouble();

        if (valorColPivote > 0) {
            double razon = valorLD / valorColPivote;
            ui->tablaSimplex->item(i, colRazon)->setText(QString::number(razon, 'f', 6));
        } else {
            ui->tablaSimplex->item(i, colRazon)->setText("-");
        }
    }
    ui->tablaSimplex->viewport()->update();
    guardarTablaActual();
    pasoActual = tablasSimplex.size() - 1;
    ui->labelPaso->setText(QString("Paso %1").arg(pasoActual));
    imprimirTablaConsola();
    actualizarFilaZ();
}


void SimplexSolver::mostrarAnalisisSensibilidad() {
    if (!ui->textEditResultado || !ui->tablaSimplex) {
        return;
    }

    using namespace Eigen;

    QString resultado;
    resultado += "ANÁLISIS DE SENSIBILIDAD\n\n";

    int filas = ui->tablaSimplex->rowCount();
    int cols = ui->tablaSimplex->columnCount();
    int colLD = cols - 2;

    QVector<QString> vbs;
    QVector<double> xB_final;
    QVector<int> vb_original_col_indices;
    int numVars = objetivo.size();
    int numHolguras_total = 0;
    int numArtificiales_total = 0;

    for (const QString& signo : signos) {
        if (signo == "<=") {
            numHolguras_total++;
        } else if (signo == ">=") {
            numHolguras_total++;
            numArtificiales_total++;
        } else if (signo == "=") {
            numArtificiales_total++;
        }
    }



    resultado += "** Variables básicas (X* = B⁻¹·b):\n";
    for (int i = 1; i < filas; ++i) {
        QString vb = ui->tablaSimplex->item(i, 1)->text().trimmed();
        double val = ui->tablaSimplex->item(i, colLD)->text().toDouble();
        vbs.append(vb);
        xB_final.append(val);
        resultado += QString("    • %1 = %2\n").arg(vb).arg(QString::number(val, 'f', 4));

        int original_col_idx = -1;
        if (vb.startsWith("x")) {
            original_col_idx = vb.mid(1).toInt() - 1;
        } else if (vb.startsWith("s")) {

            original_col_idx = numVars + (vb.mid(1).toInt() - 1);
        } else if (vb.startsWith("a")) {

            original_col_idx = numVars + numHolguras_total + (vb.mid(1).toInt() - 1);
        }

        if (original_col_idx != -1) {
            vb_original_col_indices.append(original_col_idx);
        } else {
            resultado += QString("Error: No se pudo encontrar el índice de columna original para %1\n").arg(vb);
            ui->textEditResultado->append(resultado);
            return;
        }
    }


    int n_basic = vbs.size();
    MatrixXd B(n_basic, n_basic);
    int total_vars_in_full_A = numVars + numHolguras_total + numArtificiales_total;
    QVector<QVector<double>> full_original_A_matrix(restricciones.size(), QVector<double>(total_vars_in_full_A, 0.0));

    int current_s_col_for_building_A = numVars;
    int current_a_col_for_building_A = numVars + numHolguras_total;

    for (int i = 0; i < restricciones.size(); ++i) {
        for (int j = 0; j < numVars; ++j) {
            full_original_A_matrix[i][j] = restricciones[i][j];
        }

        QString signo = signos[i];
        if (signo == "<=") {
            full_original_A_matrix[i][current_s_col_for_building_A++] = 1.0;
        } else if (signo == ">=") {
            full_original_A_matrix[i][current_s_col_for_building_A++] = -1.0;
            full_original_A_matrix[i][current_a_col_for_building_A++] = 1.0;
        } else if (signo == "=") {
            full_original_A_matrix[i][current_a_col_for_building_A++] = 1.0;
        }
    }
    for (int r = 0; r < restricciones.size(); ++r) {
        QString row_str = "  [";
        for (int c = 0; c < total_vars_in_full_A; ++c) {
            row_str += QString::number(full_original_A_matrix[r][c], 'f', 2) + " ";
        }
        row_str += "]";
    }


    for (int i = 0; i < n_basic; ++i) {
        int col_idx_in_full_original_A = vb_original_col_indices[i];
        for (int r = 0; r < n_basic; ++r) {
            B(r, i) = full_original_A_matrix[r][col_idx_in_full_original_A];
        }
    }

    resultado += "\n ** Matriz B (coeficientes de las variables básicas):\n";
    for (int i = 0; i < n_basic; ++i) {
        QString fila_str;
        for (int j = 0; j < n_basic; ++j)
            fila_str += QString::number(B(i, j), 'f', 4) + "\t";
        resultado += fila_str + "\n";
    }


    MatrixXd Binv = B.inverse();
    resultado += "\n ** Inversa de B (B⁻¹):\n";
    for (int i = 0; i < n_basic; ++i) {
        QString fila_str;
        for (int j = 0; j < n_basic; ++j)
            fila_str += QString::number(Binv(i, j), 'f', 4) + "\t";
        resultado += fila_str + "\n";
    }


    RowVectorXd cB(n_basic);
    double M = 1e6;

    for (int i = 0; i < n_basic; ++i) {
        QString vb = vbs[i];
        double original_cj_val = 0;

        if (vb.startsWith("x")) {
            original_cj_val = objetivo[vb.mid(1).toInt() - 1];
        } else if (vb.startsWith("s")) {
            original_cj_val = 0;
        } else if (vb.startsWith("a")) {
            original_cj_val = esMaximizacion ? -M : M;
        }
        cB(i) = original_cj_val;
    }
    resultado += "\n ** Vector cB (coeficientes de la base en la F.O.):\n";
    QString cB_str_output;
    for(int i=0; i<n_basic; ++i) cB_str_output += QString::number(cB(i), 'f', 4) + "\t";
    resultado += cB_str_output + "\n";


    resultado += "\n ** Precios sombra (y = c_B·B⁻¹):\n";
    RowVectorXd y_shadow_prices = cB * Binv;
    for (int i = 0; i < n_basic; ++i) {
        resultado += QString("    • Recurso %1: %2\n").arg(i + 1).arg(QString::number(y_shadow_prices(i), 'f', 4));
    }

    QVector<QString> all_tableau_var_headers;
    for (int j = 3; j < cols - 2; ++j) {
        all_tableau_var_headers.append(ui->tablaSimplex->horizontalHeaderItem(j)->text().trimmed());
    }

    resultado += "\n ** Rango permisible para coeficientes de la F.O. (Cj):\n";
    const double epsilon = 1e-8;

    for (int k_var_idx = 0; k_var_idx < all_tableau_var_headers.size(); ++k_var_idx) {
        QString var_name = all_tableau_var_headers[k_var_idx];
        double original_cj = 0;

        if (var_name.startsWith("x")) {
            original_cj = objetivo[var_name.mid(1).toInt() - 1];
        } else if (var_name.startsWith("s")) {
            original_cj = 0;
        } else if (var_name.startsWith("a")) {
            original_cj = esMaximizacion ? -M : M;
        }

        if (vbs.contains(var_name)) {
            resultado += QString("    • %1 (Variable Básica):\n").arg(var_name);
            resultado += QString("        → Cj actual = %1\n").arg(QString::number(original_cj, 'f', 2));

            double lower_bound_Ck = -std::numeric_limits<double>::infinity();
            double upper_bound_Ck = std::numeric_limits<double>::infinity();

            int row_in_tableau_for_basic_var = -1;
            for(int r = 1; r < filas; ++r) {
                if (ui->tablaSimplex->item(r, 1)->text().trimmed() == var_name) {
                    row_in_tableau_for_basic_var = r;
                    break;
                }
            }
            if (row_in_tableau_for_basic_var == -1) {
                resultado += " Error: Variable básica no encontrada en la tabla para el análisis Cj.\n";
                continue;
            }

            for (int j = 0; j < all_tableau_var_headers.size(); ++j) {
                QString nb_var_name = all_tableau_var_headers[j];

                if (!vbs.contains(nb_var_name)) {
                    double zj_cj_nb = ui->tablaSimplex->item(0, j + 3)->text().toDouble();
                    double a_bar_kj = ui->tablaSimplex->item(row_in_tableau_for_basic_var, j + 3)->text().toDouble();

                    if (std::abs(a_bar_kj) < epsilon) {
                        continue;
                    }

                    double delta_Ck_ratio = zj_cj_nb / a_bar_kj;

                    if (a_bar_kj > 0) {
                        upper_bound_Ck = std::min(upper_bound_Ck, original_cj + delta_Ck_ratio);
                    } else {
                        lower_bound_Ck = std::max(lower_bound_Ck, original_cj + delta_Ck_ratio);
                    }
                }
            }

            QString s_lb_Ck = std::isinf(lower_bound_Ck) ? "-Inf" : QString::number(lower_bound_Ck, 'f', 2);
            QString s_ub_Ck = std::isinf(upper_bound_Ck) ? "+Inf" : QString::number(upper_bound_Ck, 'f', 2);
            resultado += QString("        → Intervalo permitido: [%1, %2]\n\n")
                             .arg(s_lb_Ck)
                             .arg(s_ub_Ck);

        } else {
            resultado += QString("    • %1 (Variable No Básica):\n").arg(var_name);
            double zj_cj_current = ui->tablaSimplex->item(0, k_var_idx + 3)->text().toDouble();

            double min_cj_vnb = -std::numeric_limits<double>::infinity();
            double max_cj_vnb = original_cj + zj_cj_current;


            if (esMaximizacion && zj_cj_current < -epsilon) {
                max_cj_vnb = original_cj;
            }

            QString s_min_cj_vnb = std::isinf(min_cj_vnb) ? "-Inf" : QString::number(min_cj_vnb, 'f', 2);
            QString s_max_cj_vnb = std::isinf(max_cj_vnb) ? "+Inf" : QString::number(max_cj_vnb, 'f', 2);

            resultado += QString("        → Cj actual = %1\n").arg(QString::number(original_cj, 'f', 2));
            resultado += QString("        → Zj - Cj = %1\n").arg(QString::number(zj_cj_current, 'f', 4));
            resultado += QString("        → Intervalo permitido: [%1, %2]\n\n")
                             .arg(s_min_cj_vnb)
                             .arg(s_max_cj_vnb);
        }
    }

    resultado += "Rango permisible para lados derechos (b):\n";


    VectorXd x_basic_solution(n_basic);
    for (int i = 0; i < n_basic; ++i) {
        x_basic_solution(i) = xB_final[i];
    }

    for (int i = 0; i < n_basic; ++i) {

        double deltaMin_b = -std::numeric_limits<double>::infinity();
        double deltaMax_b = std::numeric_limits<double>::infinity();
        double bi_original = ladosDerechos[i];
        VectorXd col_of_Binv = Binv.col(i);

        for (int j = 0; j < n_basic; ++j) {
            double aji = col_of_Binv(j);
            double xj_val = x_basic_solution(j);

            if (std::abs(aji) < epsilon) {
                continue;
            }
            double ratio = -xj_val / aji;
            if (aji > 0) {
                deltaMin_b = std::max(deltaMin_b, ratio);
            } else {
                deltaMax_b = std::min(deltaMax_b, ratio);
            }
        }

        double bi_min = bi_original + deltaMin_b;
        double bi_max = bi_original + deltaMax_b;

        QString s_bi_min = std::isinf(bi_min) ? "-Inf" : QString::number(bi_min, 'f', 2);
        QString s_bi_max = std::isinf(bi_max) ? "+Inf" : QString::number(bi_max, 'f', 2);

        resultado += QString("    • b%1 (Restricción %1): [%2, %3]\n")
                         .arg(i + 1)
                         .arg(s_bi_min)
                         .arg(s_bi_max);
    }


    resultado += "\n Fórmulas utilizadas:\n";
    resultado += "    - X* = B⁻¹·b\n";
    resultado += "    - Z* = c_B·B⁻¹·b\n";
    resultado += "    - Precios sombra: y = c_B·B⁻¹\n";
    resultado += "    - Rango para Cj de VNB (Max): [-Inf, Cj_original + (Zj-Cj)_actual]\n";
    resultado += "    - Rango para Cj de VB (Max): [Cj_original + max((Zj-Cj)/a_bar_kj para a_bar_kj<0), Cj_original + min((Zj-Cj)/a_bar_kj para a_bar_kj>0)]\n";
    resultado += "    - Rango para bi: bi_original + Δ, donde Δ mantiene X* ≥ 0\n";

    ui->textEditResultado->append("\n\n-------------------------\n");
    ui->textEditResultado->append(resultado);

}

void SimplexSolver::guardarTablaActual() {
    int filas = ui->tablaSimplex->rowCount();
    int columnas = ui->tablaSimplex->columnCount();
    QVector<QString> tablaTexto;

    for (int i = 0; i < filas; ++i) {
        QString fila;
        for (int j = 0; j < columnas; ++j) {
            QTableWidgetItem* item = ui->tablaSimplex->item(i, j);
            fila += (item ? item->text() : "") + "\t";
        }
        tablaTexto.append(fila.trimmed());
    }

    if (!tablasSimplex.isEmpty() && tablaTexto == tablasSimplex.last()) {
        return;
    }

    tablasSimplex.append(tablaTexto);
}

void SimplexSolver::cargarTabla(int paso) {
    if (paso < 0 || paso >= tablasSimplex.size()) {
        return;
    }

    int filas = ui->tablaSimplex->rowCount();
    int columnas = ui->tablaSimplex->columnCount();

    QVector<QString> tablaTexto = tablasSimplex[paso];

    for (int i = 0; i < filas && i < tablaTexto.size(); ++i) {
        QStringList celdas = tablaTexto[i].split('\t');
        for (int j = 0; j < columnas && j < celdas.size(); ++j) {
            QTableWidgetItem* item = ui->tablaSimplex->item(i, j);
            if (!item) {
                item = new QTableWidgetItem();
                ui->tablaSimplex->setItem(i, j, item);
            }
            item->setText(celdas[j]);
        }
    }

    pasoActual = paso;
    actualizarLabelPaso();
    ui->tablaSimplex->viewport()->update();
}
void SimplexSolver::actualizarLabelPaso() {
    ui->labelPaso->setText(QString("Paso %1").arg(pasoActual));
}
void SimplexSolver::pasoSiguiente() {

    if (pasoActual + 1 < tablasSimplex.size()) {
        cargarTabla(pasoActual + 1);
    } else if (!solucionAlcanzada) {
        siguientePaso();
    }
}

void SimplexSolver::pasoAnterior() {
    if (pasoActual > 0) {
        cargarTabla(pasoActual - 1);
    }
}

void SimplexSolver::imprimirTablaConsola() {
    int filas = ui->tablaSimplex->rowCount();
    int columnas = ui->tablaSimplex->columnCount();

    QStringList headers;
    for (int j = 0; j < columnas; ++j) {
        headers.append(ui->tablaSimplex->horizontalHeaderItem(j)->text());
    }

    for (int i = 0; i < filas; ++i) {
        QStringList filaDatos;
        for (int j = 0; j < columnas; ++j) {
            QTableWidgetItem* item = ui->tablaSimplex->item(i, j);
            filaDatos.append(item ? item->text() : "0");
        }
    }
}

void SimplexSolver::actualizarFilaZ() {
    int cols = ui->tablaSimplex->columnCount();
    int filas = ui->tablaSimplex->rowCount();

    QVector<double> filaZ(cols, 0.0);
    double zValue = 0.0;

    for (int i = 1; i < filas; ++i) {
        QString vb = ui->tablaSimplex->item(i, 1)->text();

        if (vb.startsWith("x")) {
            bool ok;
            int idx = vb.mid(1).toInt(&ok);
            if (ok && idx >= 1 && idx <= objetivo.size()) {
                double coefVB = objetivo[idx - 1]; // c_k

                for (int j = 3; j < cols - 2; ++j) {
                    double a_kj = ui->tablaSimplex->item(i, j)->text().toDouble();
                    filaZ[j] += coefVB * a_kj;
                }

                double b_k = ui->tablaSimplex->item(i, cols - 2)->text().toDouble();
                zValue += coefVB * b_k;
            }
        }
    }

    for (int j = 3; j < cols - 2; ++j) {
        QString encabezado = ui->tablaSimplex->horizontalHeaderItem(j)->text();

        if (encabezado.startsWith("x")) {
            bool ok;
            int idx = encabezado.mid(1).toInt(&ok);
            if (ok && idx >= 1 && idx <= objetivo.size()) {
                double cj = objetivo[idx - 1];
                double zj = filaZ[j];
                double valor;
                if (esMaximizacion) {
                    valor = zj - cj;
                } else {
                    valor = cj - zj;
                }
                ui->tablaSimplex->item(0, j)->setText(QString::number(valor, 'f', 2));
                QString operacion = esMaximizacion ? "Zj - c_j" : "c_j - Zj";
            } else {
                ui->tablaSimplex->item(0, j)->setText(QString::number(filaZ[j], 'f', 2));
            }
        } else {
            ui->tablaSimplex->item(0, j)->setText(QString::number(filaZ[j], 'f', 2));
        }
    }

    ui->tablaSimplex->item(0, cols - 2)->setText(QString::number(zValue, 'f', 2));

}


void SimplexSolver::dibujarGrafico() {

    ui->stackedWidget->setCurrentWidget(ui->pagina2D);
    if (objetivo.size() != 2 || !ui->graficoCustomPlot) return;
    QCustomPlot* plot = ui->graficoCustomPlot;
    plot->clearGraphs();
    plot->clearItems();
    plot->clearPlottables();

    plot->legend->setVisible(true);
    plot->xAxis->setLabel("x1");
    plot->yAxis->setLabel("x2");
    plot->setBackground(Qt::white);
    plot->xAxis->grid()->setSubGridVisible(true);
    plot->yAxis->grid()->setSubGridVisible(true);

    int n = restricciones.size();
    QVector<QPointF> pts;

    auto esFactible = [&](double x, double y) {
        for (int i = 0; i < n; ++i) {
            double a = restricciones[i][0], b = restricciones[i][1], c = ladosDerechos[i];
            if (a * x + b * y - c > 1e-6) return false;
        }
        return x >= -1e-6 && y >= -1e-6;
    };

    for (int i = 0; i < n; ++i) {
        double a1 = restricciones[i][0], b1 = restricciones[i][1], c1 = ladosDerechos[i];
        for (int j = i + 1; j < n; ++j) {
            double a2 = restricciones[j][0], b2 = restricciones[j][1], c2 = ladosDerechos[j];
            double det = a1 * b2 - a2 * b1;
            if (fabs(det) < 1e-6) continue; // Líneas paralelas
            double xi = (c1 * b2 - c2 * b1) / det;
            double yi = (a1 * c2 - a2 * c1) / det;
            if (esFactible(xi, yi))
                pts.append(QPointF(xi, yi));
        }

        if (!qFuzzyIsNull(b1)) {
            double y = c1 / b1;
            if (esFactible(0, y)) pts.append(QPointF(0, y));
        }
        if (!qFuzzyIsNull(a1)) {
            double x = c1 / a1;
            if (esFactible(x, 0)) pts.append(QPointF(x, 0));
        }
    }

    if (esFactible(0, 0)) pts.append(QPointF(0, 0));

    QPointF centro(0, 0);
    for (auto &p : pts) centro += p;
    centro /= pts.size();
    std::sort(pts.begin(), pts.end(), [&](const QPointF &a, const QPointF &b) {
        return atan2(a.y() - centro.y(), a.x() - centro.x()) < atan2(b.y() - centro.y(), b.x() - centro.x());
    });

    if (!pts.isEmpty()) {
        QCPCurve *region = new QCPCurve(plot->xAxis, plot->yAxis);
        for (auto &p : pts) region->addData(p.x(), p.y());
        region->addData(pts[0].x(), pts[0].y());  // Cierra el polígono
        region->setBrush(QBrush(QColor(100, 200, 100, 80)));
        region->setPen(Qt::NoPen);
        region->setName("Región factible");

    }

    double maxX = 0, maxY = 0;
    for (const auto &p : pts) {
        if (p.x() > maxX) maxX = p.x();
        if (p.y() > maxY) maxY = p.y();
    }
    double rango = std::max(maxX, maxY) + 20;

    for (int i = 0; i < n; ++i) {
        QVector<double> x(2), y(2);
        double a = restricciones[i][0], b = restricciones[i][1], c = ladosDerechos[i];
        if (qFuzzyIsNull(b)) {
            double xVal = c / a;
            x[0] = x[1] = xVal;
            y[0] = 0;
            y[1] = rango;
        } else {
            x[0] = 0;
            y[0] = c / b;
            x[1] = rango;
            y[1] = (c - a * rango) / b;
        }
        QCPGraph *g = plot->addGraph();
        g->setData(x, y);
        g->setPen(QPen(QColor::fromHsv((i * 60) % 360, 200, 100), 2));
        g->setName(QString("Restricción %1").arg(i + 1));
    }

    double x1 = 0, x2 = 0;
    int colLD = ui->tablaSimplex->columnCount() - 2;
    for (int i = 1; i < ui->tablaSimplex->rowCount(); ++i) {
        QString vb = ui->tablaSimplex->item(i, 1)->text();
        double val = ui->tablaSimplex->item(i, colLD)->text().toDouble();
        if (vb == "x1") x1 = val;
        if (vb == "x2") x2 = val;
    }

    QCPGraph *punto = plot->addGraph();
    punto->setLineStyle(QCPGraph::lsNone);
    punto->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, Qt::red, 10));
    punto->addData(x1, x2);
    punto->setName("Óptimo");

    QCPItemText *texto = new QCPItemText(plot);
    texto->setPositionAlignment(Qt::AlignHCenter | Qt::AlignBottom);
    texto->position->setCoords(x1, x2 + 2.5);
    texto->setText(QString("x1 = %1\nx2 = %2").arg(x1).arg(x2));
    texto->setFont(QFont("Arial", 9));
    texto->setPen(QPen(Qt::black));
    texto->setBrush(QBrush(QColor(255, 255, 255, 180)));

    plot->xAxis->setRange(0, rango);
    plot->yAxis->setRange(0, rango);
    plot->replot();
}




void SimplexSolver::dibujarGrafico3D() {
    ui->stackedWidget->setCurrentWidget(ui->pagina3D);

    if (objetivo.size() != 3) {
        qWarning() << "The problem does not have 3 variables. 3D graph requires 3 decision variables.";
        QMessageBox::warning(nullptr, "Simplex Problem", "3D graph is only valid for problems with 3 variables (x1, x2, x3).");
        return;
    }

    QJsonArray restriccionesJson;
    for (int i = 0; i < restricciones.size(); ++i) {
        QJsonObject r;
        r["a"] = restricciones[i][0];
        r["b"] = restricciones[i][1];
        r["c"] = restricciones[i][2];
        r["rhs"] = ladosDerechos[i];
        restriccionesJson.append(r);
    }

    double x1 = 0, x2 = 0, x3 = 0;
    int colLD = ui->tablaSimplex->columnCount() - 2;

    for (int i = 1; i < ui->tablaSimplex->rowCount(); ++i) {
        QTableWidgetItem *vbItem = ui->tablaSimplex->item(i, 1);
        QTableWidgetItem *valItem = ui->tablaSimplex->item(i, colLD);

        if (vbItem && valItem) {
            QString vb = vbItem->text();
            double val = valItem->text().toDouble();
            if (vb == "x1") x1 = val;
            if (vb == "x2") x2 = val;
            if (vb == "x3") x3 = val;
        }
    }

    QJsonArray verticesJson;
    QVector<QVector3D> vertices = calcularVerticesFactibles3D();

    if (vertices.isEmpty()) {
        qWarning() << "No feasible vertices found to graph.";
        QMessageBox::warning(nullptr, "3D Graph", "No feasible vertices in the region to display.");
        return;
    }

    for (const QVector3D& v : vertices) {
        QJsonObject punto;
        punto["x"] = v.x();
        punto["y"] = v.y();
        punto["z"] = v.z();
        verticesJson.append(punto);
    }

    QJsonObject json;
    json["objetivo"] = QJsonArray{ objetivo[0], objetivo[1], objetivo[2] };
    json["restricciones"] = restriccionesJson;
    json["optimo"] = QJsonObject{ {"x", x1}, {"y", x2}, {"z", x3} };
    json["vertices"] = verticesJson;

    QJsonDocument doc(json);

    QString pythonExecutablePath;
    QString rutaScript;
    QString dataDirPath;
    QString rutaJson;

    QString appDirPath = QCoreApplication::applicationDirPath();

    QDir currentAppDir(appDirPath);
    currentAppDir.cdUp();
    currentAppDir.cdUp();

    QString projectRootPath = currentAppDir.path();

#ifdef Q_OS_WIN
    pythonExecutablePath = QDir::cleanPath(projectRootPath + "/_python_env/Scripts/python.exe");
    rutaScript = QDir::cleanPath(projectRootPath + "/grafico3D.py");
    dataDirPath = QDir::cleanPath(projectRootPath + "/data");
#else
    pythonExecutablePath = QDir::cleanPath(projectRootPath + "/_python_env/bin/python");
    rutaScript = QDir::cleanPath(projectRootPath + "/grafico3D.py");
    dataDirPath = QDir::cleanPath(projectRootPath + "/data");
#endif

    QDir().mkpath(dataDirPath);
    rutaJson = QDir::cleanPath(dataDirPath + "/datos_grafico3D.json");

    QFile archivo(rutaJson);
    if (!archivo.open(QIODevice::WriteOnly)) {
        qWarning() << "Could not write JSON file:" << rutaJson;
        QMessageBox::critical(nullptr, "File Error", "Could not write JSON data to file: " + rutaJson);
        return;
    }
    archivo.write(doc.toJson());
    archivo.close();

    qDebug() << "Calculated pythonExecutablePath:" << pythonExecutablePath;
    qDebug() << "Calculated rutaScript:" << rutaScript;
    qDebug() << "Calculated rutaJson:" << rutaJson;

    QStringList args;
    QString venvLibPath = QDir::cleanPath(projectRootPath + "/_python_env/Lib/site-packages");

    args << "-S"
         << "-E"
         << "-c"
         << QString("import sys; sys.path.insert(0, '%1'); import runpy; runpy.run_path('%2', run_name='__main__', init_globals=globals())")
                .arg(venvLibPath)
                .arg(rutaScript)
         << rutaJson;

    QProcess *proceso = new QProcess(this);

    QProcessEnvironment env = QProcessEnvironment::systemEnvironment();
#ifdef Q_OS_WIN
    QString currentPath = env.value("PATH");
    QString pythonDLLPath = QDir::cleanPath(projectRootPath + "/_python_env");
    QString pythonScriptsPath = QDir::cleanPath(projectRootPath + "/_python_env/Scripts");
    env.insert("PATH", pythonDLLPath + ";" + pythonScriptsPath + ";" + currentPath);
#else
    QString currentPath = env.value("PATH");
    QString pythonBinPath = QDir::cleanPath(projectRootPath + "/_python_env/bin");
    env.insert("PATH", pythonBinPath + ":" + currentPath);
#endif
    qDebug() << "Process environment PATH set to:" << env.value("PATH");
    proceso->setProcessEnvironment(env);

    connect(proceso, &QProcess::readyReadStandardOutput, [proceso]() {
        qDebug() << "stdout (Python):" << proceso->readAllStandardOutput();
    });
    connect(proceso, &QProcess::readyReadStandardError, [proceso]() {
        qDebug() << "stderr (Python):" << proceso->readAllStandardError();
    });
    connect(proceso, QOverload<int, QProcess::ExitStatus>::of(&QProcess::finished),
            [proceso](int exitCode, QProcess::ExitStatus status) {
                qDebug() << "Python process finished with code:" << exitCode << "status:" << status;
                if (exitCode != 0) {
                    QMessageBox::critical(nullptr, "Python Script Error", "The 3D graph script encountered an error. Check application logs for details.");
                }
                proceso->deleteLater();
            });

    proceso->start(pythonExecutablePath, args);
    if (!proceso->waitForStarted(5000)) {
        qCritical() << "Failed to start Python process at:" << pythonExecutablePath;
        QMessageBox::critical(nullptr, "Process Error", "Failed to start Python script. Check if Python environment is correctly packaged and path is correct.");
    }
}


QVector<QVector3D> SimplexSolver::calcularVerticesFactibles3D() {
    QVector<QVector3D> vertices;
    int n = restricciones.size();

    auto cumpleRestricciones = [&](const Eigen::Vector3d& punto) -> bool {
        for (int m = 0; m < n; ++m) {
            double a = restricciones[m][0];
            double b = restricciones[m][1];
            double c = restricciones[m][2];
            double rhs = ladosDerechos[m];
            if (a * punto(0) + b * punto(1) + c * punto(2) - rhs > 1e-2)
                return false;
        }
        return punto(0) >= -1e-5 && punto(1) >= -1e-5 && punto(2) >= -1e-5;
    };

    auto agregarSiNoDuplicado = [&](const Eigen::Vector3d& punto) {
        for (const QVector3D& v : vertices) {
            if ((std::abs(v.x() - punto(0)) < 1e-3) &&
                (std::abs(v.y() - punto(1)) < 1e-3) &&
                (std::abs(v.z() - punto(2)) < 1e-3)) {
                return;
            }
        }
        vertices.append(QVector3D(punto(0), punto(1), punto(2)));
    };

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            for (int k = j + 1; k < n; ++k) {
                Eigen::Matrix3d A;
                A << restricciones[i][0], restricciones[i][1], restricciones[i][2],
                    restricciones[j][0], restricciones[j][1], restricciones[j][2],
                    restricciones[k][0], restricciones[k][1], restricciones[k][2];
                Eigen::Vector3d b(ladosDerechos[i], ladosDerechos[j], ladosDerechos[k]);

                if (std::abs(A.determinant()) < 1e-8) continue;

                Eigen::Vector3d punto = A.colPivHouseholderQr().solve(b);
                if (cumpleRestricciones(punto))
                    agregarSiNoDuplicado(punto);
            }
        }
    }

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            for (int varFija = 0; varFija < 3; ++varFija) {
                Eigen::Matrix3d A;
                Eigen::Vector3d b;

                A.row(0) = Eigen::Vector3d::Zero();
                A.row(0)[varFija] = 1;
                b(0) = 0;

                A.row(1) = Eigen::Vector3d(restricciones[i][0], restricciones[i][1], restricciones[i][2]);
                A.row(2) = Eigen::Vector3d(restricciones[j][0], restricciones[j][1], restricciones[j][2]);

                b(1) = ladosDerechos[i];
                b(2) = ladosDerechos[j];

                if (std::abs(A.determinant()) < 1e-8) continue;

                Eigen::Vector3d punto = A.colPivHouseholderQr().solve(b);
                if (cumpleRestricciones(punto))
                    agregarSiNoDuplicado(punto);
            }
        }
    }

    for (int i = 0; i < n; ++i) {
        for (int var1 = 0; var1 < 3; ++var1) {
            for (int var2 = var1 + 1; var2 < 3; ++var2) {
                Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
                Eigen::Vector3d b = Eigen::Vector3d::Zero();

                A(var1, var1) = 1;
                A(var2, var2) = 1;
                A(3 - var1 - var2, 0) = restricciones[i][0];
                A(3 - var1 - var2, 1) = restricciones[i][1];
                A(3 - var1 - var2, 2) = restricciones[i][2];
                b(3 - var1 - var2) = ladosDerechos[i];

                if (std::abs(A.determinant()) < 1e-8) continue;

                Eigen::Vector3d punto = A.colPivHouseholderQr().solve(b);
                if (cumpleRestricciones(punto))
                    agregarSiNoDuplicado(punto);
            }
        }
    }

    return vertices;
}



void SimplexSolver::reiniciarEstado() {
    tablasSimplex.clear();
    pasoActual = 0;
    solucionAlcanzada = false;
    columnaPivote = -1;


    if (ui->textEditResultado) {
        ui->textEditResultado->clear();
    }


    if (ui->graficoCustomPlot) {
        ui->graficoCustomPlot->clearGraphs();
        ui->graficoCustomPlot->clearItems();
        ui->graficoCustomPlot->clearPlottables();
        ui->graficoCustomPlot->replot();
    }

}

