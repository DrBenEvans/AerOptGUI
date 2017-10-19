#include "windowsizing.h"
#include <QApplication>
#include <QDebug>
#include <QStyle>

void centerAndResizeWindow(QWidget* window, float width, float height) {
    // get the dimension available on this screen
    QSize availableSize = qApp->desktop()->availableGeometry().size();
    width *= availableSize.width();
    height *= availableSize.height();
    qDebug() << "Available window dimensions " << width << "x" << height;
    qDebug() << "Computed window dimensions " << width << "x" << height;
    QSize newSize( width, height );

    window->setGeometry(
        QStyle::alignedRect(
            Qt::LeftToRight,
            Qt::AlignCenter,
            newSize,
            qApp->desktop()->availableGeometry()
        )
    );
}
