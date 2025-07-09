#include "mainwindow.h"
#include <QApplication>
#include <QFile>
#include <QDebug>
#include <QFont>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    QFont font("Segoe UI", 10);
    a.setFont(font);

    QFile styleFile(":/styles.qss");
    if (styleFile.open(QFile::ReadOnly | QFile::Text)) {
        QString style = styleFile.readAll();
        a.setStyleSheet(style);
        styleFile.close();
        qDebug() << "Stylesheet loaded successfully.";
    } else {
        qDebug() << " Could not open stylesheet file. Make sure it's in resources.qrc and path is correct: :/styles/styles.qss";
    }

    MainWindow w;
    w.show();
    return a.exec();
}
