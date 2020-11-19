/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Charts module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#include "ChartView.h"
#include <QtGui/QMouseEvent>
#include <QApplication>

ChartView::ChartView(QWidget *parent) :
    QChartView(parent)
{
    setDragMode(QGraphicsView::NoDrag);
    this->setMouseTracking(true);

}

void ChartView::wheelEvent(QWheelEvent* event) {
    if(event->delta() > 0) {
        chart()->zoomIn();
    } else {
        chart()->zoomOut();
    }
}

//Credit: https://stackoverflow.com/questions/46805186/qt-chart-move-view-with-pressed-middle-mouse-button
void ChartView::mousePressEvent(QMouseEvent *event)
{
    if (event->button() == Qt::LeftButton)
    {
        QApplication::setOverrideCursor(QCursor(Qt::SizeAllCursor));
        m_lastMousePos = event->pos();
        event->accept();
    }

    QChartView::mousePressEvent(event);
}

//Credit: https://stackoverflow.com/questions/46805186/qt-chart-move-view-with-pressed-middle-mouse-button
void ChartView::mouseMoveEvent(QMouseEvent *event)
{
    // pan the chart with a middle mouse drag
    if (event->buttons() & Qt::LeftButton)
    {
        auto dPos = event->pos() - m_lastMousePos;
        chart()->scroll(-dPos.x(), dPos.y());

        m_lastMousePos = event->pos();
        event->accept();

        QApplication::restoreOverrideCursor();
    }

    QChartView::mouseMoveEvent(event);
}
