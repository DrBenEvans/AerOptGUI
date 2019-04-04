#include "MeshView.h"
#include "BoundaryPointView.h"
#include <QtWidgets>

MeshView::MeshView(ViewScaler* scale, QGraphicsItem* parent) :
    QGraphicsObject(parent),
    mMeshModel(nullptr),
    mScale(scale),
    mBoundaryPointModel(nullptr)
{
    this->color = QColor(0,255,0);
    setZValue(-10);
    setPos(0,0);

    // set initial bounding box
    mXmax =  -std::numeric_limits<float>::infinity();
    mYmax =  -std::numeric_limits<float>::infinity();
    mXmin =  std::numeric_limits<float>::infinity();
    mYmin =  std::numeric_limits<float>::infinity();
}

void MeshView::setMeshModel(MeshModel* model) {
    mMeshModel = model;
    connect(mMeshModel, &MeshModel::meshChanged, this, &MeshView::meshChanged);
    meshChanged();
}

void MeshView::setBoundaryPointModel(BoundaryPointModel* model) {
    mBoundaryPointModel = model;
    connect(mBoundaryPointModel, &BoundaryPointModel::boundaryPointsReset, this, &MeshView::boundaryPointsReset);
    boundaryPointsReset();
}

void MeshView::meshChanged() {
    prepareGeometryChange();
    calcBoundingBox();
}

void MeshView::calcBoundingBox() {
    const Mesh* mesh = mMeshModel->currentMesh();
    if(!mesh) return;

    float x, y;

    for (auto& p : mesh->getMeshPoints())
    {
        x = mScale->w(p.first);
        y = mScale->h(p.second);

        if(x>mXmax)
            mXmax = x;
        if(x<mXmin)
            mXmin = x;
        if(y>mYmax)
            mYmax = y;
        if(y<mYmin)
            mYmin = y;

    }

    int margin=10;

    mBoundingBox = QRectF(mXmin-margin, mYmin-margin, (mXmax-mXmin)+2*margin, (mYmax-mYmin)+2*margin);

    QGraphicsScene * qscene = scene();
    if( qscene ){
        qscene->setSceneRect(mBoundingBox);
    }
}

void MeshView::boundaryPointsReset()
{
    foreach(auto& item, this->childItems()) {
        delete item;
    }

    // Draw control node objects, if mesh is set
    if(mBoundaryPointModel) {
        for (int index=0; index < mBoundaryPointModel->count(); index++)
        {
            BoundaryPointView* bpView = new BoundaryPointView(mBoundaryPointModel, index, mScale, this);
        }
    }
    update();
}

QRectF MeshView::boundingRect() const
{
    return mBoundingBox;
}

QPainterPath MeshView::shape() const
{
    QPainterPath path;
    path.addRect(boundingRect());
    return path;
}

void MeshView::paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
    Q_UNUSED(widget);
    const Mesh* mesh = mMeshModel->currentMesh();
    if(!mesh) return;

    painter->setRenderHint(QPainter::Antialiasing, false);
    painter->setBrush(QBrush(QColor(255, 255, 255, 0)));
    painter->setPen(QPen(Qt::black, getBrushSize(), Qt::SolidLine));
    QRect rect;
    rect.setWidth( 7 );
    rect.setHeight( 7 );

    const auto& bConnects = mesh->getBConnects();
    const auto& meshpoints = mesh->getMeshPoints();
    const auto& meshconnects = mesh->getMeshConnectivities();
    const auto& resultsdata = mesh->getMeshData();

    //This scope is for rendering the mesh when available
    painter->setPen( QPen(Qt::lightGray, getBrushSize(), Qt::SolidLine) );
    painter->setRenderHint(QPainter::Antialiasing, true);
    if (meshconnects.size() > 0)
    {
        const uint r = 4;//pressure, 0 = rho, 1 = u, 2 = v, 3 = energy

        //Determine and set scalar range
        float min = std::numeric_limits<float>::infinity();
        float max = -std::numeric_limits<float>::infinity();
        for (const auto& s : resultsdata)
        {
            const float& v = std::get<r>(s);
            if (v < min) min = v;
            if (v > max) max = v;
        }
        mColorMapper.setMin(min);
        mColorMapper.setMax(max);

        //Render triangulation
        QColor color;

        for (auto& tuple : meshconnects)
        {
            const uint& i = std::get<0>(tuple);
            const uint& j = std::get<1>(tuple);
            const uint& k = std::get<2>(tuple);

            const std::pair<float,float>& p1 = meshpoints.at(i-1);
            const std::pair<float,float>& p2 = meshpoints.at(j-1);
            const std::pair<float,float>& p3 = meshpoints.at(k-1);

            if (resultsdata.size() == meshpoints.size())
            {
                const float& pressure1 = std::get<r>( resultsdata.at(i-1) );//Pressure
                const float& pressure2 = std::get<r>( resultsdata.at(j-1) );//Pressure
                const float& pressure3 = std::get<r>( resultsdata.at(k-1) );//Pressure
                const float pressure = (pressure1 + pressure2 + pressure3) * 0.333;

                color = mColorMapper.color(pressure);

                QPainterPath path;
                path.moveTo( qreal(mScale->w( p1.first )), qreal(mScale->h( p1.second )) );
                path.lineTo( qreal(mScale->w( p2.first )), qreal(mScale->h( p2.second )) );
                path.lineTo( qreal(mScale->w( p3.first )), qreal(mScale->h( p3.second )) );
                path.lineTo( qreal(mScale->w( p1.first )), qreal(mScale->h( p1.second )) );
                painter->fillPath(path, QBrush(color));
            }

            //Draw triangle boundary here
            painter->drawLine(
                        mScale->w( p1.first ),
                        mScale->h( p1.second ),
                        mScale->w( p2.first ),
                        mScale->h( p2.second )
                        );

            painter->drawLine(
                        mScale->w( p2.first ),
                        mScale->h( p2.second ),
                        mScale->w( p3.first ),
                        mScale->h( p3.second )
                        );

            painter->drawLine(
                        mScale->w( p3.first ),
                        mScale->h( p3.second ),
                        mScale->w( p1.first ),
                        mScale->h( p1.second )
                        );
        }
    }
}

uint MeshView::getBrushSize() {
    return 0;
}
