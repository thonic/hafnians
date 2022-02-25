pipeline {
    agent any
    stages {
        stage('Run hafnian.py') {
            steps {
                sh '''
                    pwd
                    whoami
                    conda env list
                '''
            }
        }
    }
}