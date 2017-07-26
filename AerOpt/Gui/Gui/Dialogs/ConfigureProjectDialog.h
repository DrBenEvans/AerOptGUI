#ifndef CONFIGUREPROJECT_H
#define CONFIGUREPROJECT_H

#include <QDialog>

namespace Ui {
class ConfigureProject;
}

class ConfigureProject : public QDialog
{
    Q_OBJECT

public:
    explicit ConfigureProject(QWidget *parent = 0);
    ~ConfigureProject();

private:
    Ui::ConfigureProject *ui;
};

#endif // CONFIGUREPROJECT_H
